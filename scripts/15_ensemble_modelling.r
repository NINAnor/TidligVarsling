# script to fit GAM and boosted regression tree models for plants and insects
# test predictive ability and generate weights for ensemble models
# Jenny Hansen
# 29 September 2025
# updated 02 October 2025 to try with 1000 bg pts
# updated 08 October 2025 with new data & adding MaxEnt
# updated 05 Nov 2025 to include new plant data
# will also use all variables instead of select variables

# NB! This script takes 3+ hours to run

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(mgcv)
library(gbm)
library(pROC)
library(PresenceAbsence)
library(blockCV)
library(terra)
library(tibble)
library(randomForest)
library(maxnet)
library(parallel)


# Import data -------------------------------------------------------------

# GAM, BRT, and RF data (script 13)
insect_sf <- st_read("vector/insect_model_sf.geojson")  %>% 
  mutate(group = "insect", present = as.integer(species_richness > 0))
plant_sf  <- st_read("vector/plant_model_sf.geojson") %>% 
  mutate(group = "plant", present = as.integer(species_richness > 0))

# maxent data (script 13)
insect_me <- st_read("vector/insect_maxent_data.geojson")
plant_me <- st_read("vector/plant_maxent_data.geojson")

insect_me_df <- insect_me %>% st_drop_geometry() %>% tidyr::drop_na()
plant_me_df <- plant_me %>% st_drop_geometry() %>% tidyr::drop_na()

# predictors from script 10
pred <- rast("raster/complete_prediction_stack.tif")

# Set up spatial folds ----------------------------------------------------

set.seed(42)

# determine spatial autocorrelation range

# using point data (insects)
cv_range_points <- cv_spatial_autocor(x = insect_sf, column = "present", 
                                      plot = TRUE)
cv_range_points$range # 2326.377

# using point data (plants)
cv_range_points <- cv_spatial_autocor(x = plant_sf, column = "present", 
                                      plot = TRUE)
cv_range_points$range # 24681.96

# set block size
block_i = 2327
block_p <- 24682

folds_insects <- cv_spatial(x = insect_sf, "present", k = 5, size = block_i,
                            selection = "random", iteration = 200)
folds_plants <- cv_spatial(x = plant_sf, "present", k = 5, size = block_p,
                           selection = "random", iteration = 200)


# GAM hurdle- insects -----------------------------------------------------

# set up data
preds_insects <- names(insect_sf %>% 
                         select(distance_to_public_road:percentage_agriculture) %>% 
                         st_drop_geometry())
dat_insects <- insect_sf %>% st_drop_geometry()
dat_insects$present <- as.integer(dat_insects$species_richness > 0)

# out-of-fold storage
oof_prob_insects <- rep(NA_real_, nrow(dat_insects))  # probability of >0
oof_mu_insects   <- rep(NA_real_, nrow(dat_insects))  # conditional mean richness


# 2-part hurdle model 
for (i in seq_along(folds_insects$folds_list)) {
  idx_train <- folds_insects$folds_list[[i]][[1]]  
  idx_test  <- folds_insects$folds_list[[i]][[2]]  
  
  train <- dat_insects[idx_train, ]
  test  <- dat_insects[idx_test, ]
  
  # occurrence model (binomial gam)
  f_bin <- as.formula(paste("present ~",
                            paste(sprintf("s(%s,k=4,bs='cs')", preds_insects), 
                                  collapse=" + ")))
  gam_bin <- mgcv::gam(f_bin, family=binomial(), data=train, method="REML")
  
  oof_prob_insects[idx_test] <- predict(gam_bin, newdata=test, type="response")
  
  # conditional richness model (Gaussian gam)
  train_pos <- train[train$present == 1, ]
  if (nrow(train_pos) >= 10) {
    f_mu <- as.formula(paste("species_richness ~",
                             paste(sprintf("s(%s,k=4,bs='cs')", preds_insects), 
                                   collapse=" + ")))
    gam_mu <- mgcv::gam(f_mu, family=gaussian(), data=train_pos, method="REML")
    oof_mu_insects[idx_test] <- predict(gam_mu, newdata=test, type="response")
  } else {
    oof_mu_insects[idx_test] <- NA_real_
  }
}

obs_insects <- dat_insects$species_richness
expected_richness_insects <- oof_prob_insects * oof_mu_insects

#metrics on conditional part
rmse_exp_insects <- sqrt(mean((expected_richness_insects - obs_insects)^2, 
                              na.rm=TRUE))
r2_exp_insects   <- 1 - sum((obs_insects - expected_richness_insects)^2, 
                            na.rm=TRUE) /
  sum((obs_insects - mean(obs_insects, na.rm=TRUE))^2, na.rm=TRUE)

# metrics on occurrence part
y_insects <- dat_insects$present
prob_insects <- oof_prob_insects

auc_occ_insects <- pROC::auc(y_insects, prob_insects)

thresh <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:length(y_insects), observed=y_insects, predicted=prob_insects),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm <- PresenceAbsence::cmx(
  data.frame(id=1:length(y_insects), observed=y_insects, predicted=prob_insects),
  threshold=thresh)

tss_occ_insects <- PresenceAbsence::sensitivity(cm) + 
  PresenceAbsence::specificity(cm) - 1

# GAM hurdle- plants ------------------------------------------------------

# set up data
preds_plants <- names(plant_sf %>% 
                        select(distance_to_public_road:percentage_agriculture) %>% 
                        st_drop_geometry())
dat_plants <- plant_sf %>% st_drop_geometry()
dat_plants$present <- as.integer(dat_plants$species_richness > 0)

# out-of-fold storage
oof_prob_plants <- rep(NA_real_, nrow(dat_plants))  
oof_mu_plants   <- rep(NA_real_, nrow(dat_plants))  


# 2-part hurdle model 
for (i in seq_along(folds_plants$folds_list)) {
  idx_train <- folds_plants$folds_list[[i]][[1]] 
  idx_test  <- folds_plants$folds_list[[i]][[2]]  
  
  train <- dat_plants[idx_train, ]
  test  <- dat_plants[idx_test, ]

  f_bin <- as.formula(paste("present ~",
                            paste(sprintf("s(%s,k=4,bs='cs')", preds_plants), 
                                  collapse=" + ")))
  gam_bin <- mgcv::gam(f_bin, family=binomial(), data=train, method="REML")
  
  oof_prob_plants[idx_test] <- predict(gam_bin, newdata=test, type="response")
  
  train_pos <- train[train$present == 1, ]
  if (nrow(train_pos) >= 10) {
    f_mu <- as.formula(paste("species_richness ~",
                             paste(sprintf("s(%s,k=4,bs='cs')", preds_plants), 
                                   collapse=" + ")))
    gam_mu <- mgcv::gam(f_mu, family=gaussian(), data=train_pos, method="REML")
    oof_mu_plants[idx_test] <- predict(gam_mu, newdata=test, type="response")
  } else {
    oof_mu_plants[idx_test] <- NA_real_
  }
}

obs_plants <- dat_plants$species_richness
expected_richness_plants <- oof_prob_plants * oof_mu_plants

rmse_exp_plants <- sqrt(mean((expected_richness_plants - obs_plants)^2, 
                             na.rm=TRUE))
r2_exp_plants   <- 1 - sum((obs_plants - expected_richness_plants)^2, 
                           na.rm=TRUE) /
  sum((obs_plants - mean(obs_plants, na.rm=TRUE))^2, na.rm=TRUE)

y_plants <- dat_plants$present
prob_plants <- oof_prob_plants

auc_occ_plants <- pROC::auc(y_plants, prob_plants)

thresh <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:length(y_plants), observed=y_plants, predicted=prob_plants),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm <- PresenceAbsence::cmx(
  data.frame(id=1:length(y_plants), observed=y_plants, predicted=prob_plants),
  threshold=thresh)

tss_occ_plants <- PresenceAbsence::sensitivity(cm) + 
  PresenceAbsence::specificity(cm) - 1

# BRT hurdle -insects -----------------------------------------------------

# setup insects
oof_prob_brt_insects <- rep(NA_real_, nrow(dat_insects))  # probability of >0
oof_mu_brt_insects   <- rep(NA_real_, nrow(dat_insects))  # conditional mean richness


# 2-part BRT hurdle model
for (i in seq_along(folds_insects$folds_list)) {
  idx_train <- folds_insects$folds_list[[i]][[1]]
  idx_test  <- folds_insects$folds_list[[i]][[2]]
  
  train <- dat_insects[idx_train, ]
  test  <- dat_insects[idx_test, ]
  
  # occurrence model (bernoulli) 
  df_occ <- data.frame(present = train$present, train[, preds_insects, 
                                                      drop=FALSE])
  
  brt_bin <- gbm::gbm(
    formula = present ~ .,
    data = df_occ,
    distribution = "bernoulli",
    n.trees = 5000, interaction.depth = 3,
    shrinkage = 0.01, bag.fraction = 0.6,
    cv.folds = 5, verbose = FALSE
  )
  best_iter_bin <- gbm::gbm.perf(brt_bin, method="cv", plot.it=FALSE)
  
  oof_prob_brt_insects[idx_test] <- predict(
    brt_bin,
    newdata = test[, preds_insects, drop=FALSE],
    n.trees = best_iter_bin,
    type = "response"
  )
  
  # conditional richness model (Gaussian)
  train_pos <- train[train$present == 1, ]
  if (nrow(train_pos) >= 10) {
    df_mu <- data.frame(species_richness = train_pos$species_richness,
                        train_pos[, preds_insects, drop=FALSE])
    
    brt_mu <- gbm::gbm(
      formula = species_richness ~ .,
      data = df_mu,
      distribution = "gaussian",
      n.trees = 5000, interaction.depth = 3,
      shrinkage = 0.01, bag.fraction = 0.6,
      cv.folds = 5, verbose = FALSE
    )
    best_iter_mu <- gbm::gbm.perf(brt_mu, method="cv", plot.it=FALSE)
    
    oof_mu_brt_insects[idx_test] <- predict(
      brt_mu,
      newdata = test[, preds_insects, drop=FALSE],
      n.trees = best_iter_mu,
      type = "response"
    )
  } else {
    oof_mu_brt_insects[idx_test] <- NA_real_
  }
}

# performance metrics
expected_richness_brt_insects <- oof_prob_brt_insects * oof_mu_brt_insects
rmse_brt_insects <- sqrt(mean((expected_richness_brt_insects - 
                                 dat_insects$species_richness)^2, 
                              na.rm=TRUE))
r2_brt_insects <- 1 - sum((dat_insects$species_richness - 
                             expected_richness_brt_insects)^2, 
                          na.rm=TRUE) /
  sum((dat_insects$species_richness - mean(dat_insects$species_richness))^2, 
      na.rm=TRUE)

y_occ_brt_insects <- dat_insects$present
prob_occ_brt_insects <- oof_prob_brt_insects

# AUC
auc_occ_brt_insects <- pROC::auc(y_occ_brt_insects, 
                                 prob_occ_brt_insects)

thresh_brt_insects <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:length(y_occ_brt_insects),
             observed=y_occ_brt_insects,
             predicted=prob_occ_brt_insects),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm_brt_insects <- PresenceAbsence::cmx(
  data.frame(id=1:length(y_occ_brt_insects),
             observed=y_occ_brt_insects,
             predicted=prob_occ_brt_insects),
  threshold=thresh_brt_insects)

tss_occ_brt_insects <- PresenceAbsence::sensitivity(cm_brt_insects) +
  PresenceAbsence::specificity(cm_brt_insects) - 1

# BRT hurdle- plants ------------------------------------------------------

# setup
oof_prob_brt_plants <- rep(NA_real_, nrow(dat_plants))  
oof_mu_brt_plants   <- rep(NA_real_, nrow(dat_plants))  


# 2-part BRT hurdle model
for (i in seq_along(folds_plants$folds_list)) {
  idx_train <- folds_plants$folds_list[[i]][[1]]
  idx_test  <- folds_plants$folds_list[[i]][[2]]
  
  train <- dat_plants[idx_train, ]
  test  <- dat_plants[idx_test, ]

  df_occ <- data.frame(present = train$present, train[, preds_plants, 
                                                      drop=FALSE])

  brt_bin <- gbm::gbm(
    formula = present ~ .,
    data = df_occ,
    distribution = "bernoulli",
    n.trees = 5000, interaction.depth = 2,
    shrinkage = 0.01, bag.fraction = 0.6,
    n.minobsinnode = 3,
    cv.folds = 5, verbose = FALSE
  )
  best_iter_bin <- gbm::gbm.perf(brt_bin, method="cv", plot.it=FALSE)
  
  oof_prob_brt_plants[idx_test] <- predict(
    brt_bin,
    newdata = test[, preds_plants, drop=FALSE],
    n.trees = best_iter_bin,
    type = "response"
  )
  
  train_pos <- train[train$present == 1, ]
  if (nrow(train_pos) >= 10) {
    df_mu <- data.frame(species_richness = train_pos$species_richness,
                        train_pos[, preds_plants, drop=FALSE])
    
    brt_mu <- gbm::gbm(
      formula = species_richness ~ .,
      data = df_mu,
      distribution = "gaussian",
      n.trees = 5000, interaction.depth = 2,
      shrinkage = 0.01,
      bag.fraction = 0.6,
      n.minobsinnode = 3,
      cv.folds = 5,
      verbose = FALSE
    )
    best_iter_mu <- gbm::gbm.perf(brt_mu, method="cv", plot.it=FALSE)
    
    oof_mu_brt_plants[idx_test] <- predict(
      brt_mu,
      newdata = test[, preds_plants, drop=FALSE],
      n.trees = best_iter_mu,
      type = "response"
    )
  } else {
    oof_mu_brt_plants[idx_test] <- NA_real_
  }
}

# performance metrics
expected_richness_brt_plants <- oof_prob_brt_plants * oof_mu_brt_plants
rmse_brt_plants <- sqrt(mean((expected_richness_brt_plants - 
                                dat_plants$species_richness)^2, 
                             na.rm=TRUE))
r2_brt_plants   <- 1 - sum((dat_plants$species_richness - 
                              expected_richness_brt_plants)^2, 
                           na.rm=TRUE) /
  sum((dat_plants$species_richness - mean(dat_plants$species_richness))^2, 
      na.rm=TRUE)

y_occ_brt_plants <- dat_plants$present
prob_occ_brt_plants <- oof_prob_brt_plants


auc_occ_brt_plants <- pROC::auc(y_occ_brt_plants, 
                                prob_occ_brt_plants)

thresh_brt_plants <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:length(y_occ_brt_plants),
             observed=y_occ_brt_plants,
             predicted=prob_occ_brt_plants),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm_brt_plants <- PresenceAbsence::cmx(
  data.frame(id=1:length(y_occ_brt_plants),
             observed=y_occ_brt_plants,
             predicted=prob_occ_brt_plants),
  threshold=thresh_brt_plants)

tss_occ_brt_plants <- PresenceAbsence::sensitivity(cm_brt_plants) +
  PresenceAbsence::specificity(cm_brt_plants) - 1

# RF hurdle- insects ------------------------------------------------------

oof_prob_rf_insects <- rep(NA_real_, nrow(dat_insects))  #
oof_mu_rf_insects   <- rep(NA_real_, nrow(dat_insects))  

for (i in seq_along(folds_insects$folds_list)) {
  idx_train <- folds_insects$folds_list[[i]][[1]]
  idx_test  <- folds_insects$folds_list[[i]][[2]]
  
  train <- dat_insects[idx_train, ]
  test  <- dat_insects[idx_test, ]
  
  # occurrence model (classification RF)
  rf_bin <- randomForest(
    x = train[, preds_insects, drop=FALSE],
    y = as.factor(train$present),
    ntree = 1000, mtry = floor(sqrt(length(preds_insects)))
  )
  
  oof_prob_rf_insects[idx_test] <- predict(
    rf_bin, newdata = test[, preds_insects, drop=FALSE],
    type = "prob")[, "1"]
  
  # conditional richness model (regression RF)
  train_pos <- train[train$present == 1, ]
  if (nrow(train_pos) >= 10) {
    rf_mu <- randomForest(
      x = train_pos[, preds_insects, drop=FALSE],
      y = train_pos$species_richness,
      ntree = 1000, mtry = floor(sqrt(length(preds_insects)))
    )
    oof_mu_rf_insects[idx_test] <- predict(
      rf_mu, newdata = test[, preds_insects, drop=FALSE]
    )
  } else {
    oof_mu_rf_insects[idx_test] <- NA_real_
  }
}

# performance metrics
expected_richness_rf_insects <- oof_prob_rf_insects * oof_mu_rf_insects
rmse_rf_insects <- sqrt(mean((expected_richness_rf_insects - 
                                dat_insects$species_richness)^2, na.rm=TRUE))
r2_rf_insects   <- 1 - sum((dat_insects$species_richness - 
                              expected_richness_rf_insects)^2, na.rm=TRUE) /
  sum((dat_insects$species_richness - mean(dat_insects$species_richness))^2, 
      na.rm=TRUE)

auc_rf_insects <- pROC::auc(dat_insects$present, oof_prob_rf_insects)

thresh_rf_insects <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:length(dat_insects$present), 
             observed=dat_insects$present, 
             predicted=oof_prob_rf_insects),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm_rf_insects <- PresenceAbsence::cmx(
  data.frame(id=1:length(dat_insects$present), 
             observed=dat_insects$present, 
             predicted=oof_prob_rf_insects),
  threshold=thresh_rf_insects)

tss_rf_insects <- PresenceAbsence::sensitivity(cm_rf_insects) +
  PresenceAbsence::specificity(cm_rf_insects) - 1

# RF hurdle-plants --------------------------------------------------------

oof_prob_rf_plants <- rep(NA_real_, nrow(dat_plants))
oof_mu_rf_plants   <- rep(NA_real_, nrow(dat_plants))

for (i in seq_along(folds_plants$folds_list)) {
  idx_train <- folds_plants$folds_list[[i]][[1]]
  idx_test  <- folds_plants$folds_list[[i]][[2]]
  
  train <- dat_plants[idx_train, ]
  test  <- dat_plants[idx_test, ]

  rf_bin <- randomForest(
    x = train[, preds_plants, drop=FALSE],
    y = as.factor(train$present),
    ntree = 1000, mtry = floor(sqrt(length(preds_plants)))
  )
  
  oof_prob_rf_plants[idx_test] <- predict(
    rf_bin, newdata = test[, preds_plants, drop=FALSE],
    type = "prob")[, "1"]
  
  train_pos <- train[train$present == 1, ]
  if (nrow(train_pos) >= 10) {
    rf_mu <- randomForest(
      x = train_pos[, preds_plants, drop=FALSE],
      y = train_pos$species_richness,
      ntree = 1000, mtry = floor(sqrt(length(preds_plants)))
    )
    oof_mu_rf_plants[idx_test] <- predict(
      rf_mu, newdata = test[, preds_plants, drop=FALSE]
    )
  } else {
    oof_mu_rf_plants[idx_test] <- NA_real_
  }
}

# performance metrics
expected_richness_rf_plants <- oof_prob_rf_plants * oof_mu_rf_plants
rmse_rf_plants <- sqrt(mean((expected_richness_rf_plants - 
                               dat_plants$species_richness)^2, na.rm=TRUE))
r2_rf_plants   <- 1 - sum((dat_plants$species_richness - 
                             expected_richness_rf_plants)^2, na.rm=TRUE) /
  sum((dat_plants$species_richness - mean(dat_plants$species_richness))^2, 
      na.rm=TRUE)

auc_rf_plants <- pROC::auc(dat_plants$present, oof_prob_rf_plants)

thresh_rf_plants <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:length(dat_plants$present), 
             observed=dat_plants$present, 
             predicted=oof_prob_rf_plants),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm_rf_plants <- PresenceAbsence::cmx(
  data.frame(id=1:length(dat_plants$present), 
             observed=dat_plants$present, 
             predicted=oof_prob_rf_plants),
  threshold=thresh_rf_plants)

tss_rf_plants <- PresenceAbsence::sensitivity(cm_rf_plants) +
  PresenceAbsence::specificity(cm_rf_plants) - 1

# Run MaxEnt- insects -----------------------------------------------------

# prep data
X_insects <- insect_me_df[, -1]   
y_insects <- insect_me_df$presence

# initialize storage for out-of-fold predictions
oof_prob_insects_me <- rep(NA_real_, nrow(insect_me_df))

for (i in seq_along(folds_insects$folds_list)) {
  idx_train <- folds_insects$folds_list[[i]][[1]]
  idx_test  <- folds_insects$folds_list[[i]][[2]]
  
  m_insect_fold <- maxnet(
    p = y_insects[idx_train],
    data = X_insects[idx_train, ],
    f = maxnet.formula(p = y_insects[idx_train],
                       data = X_insects[idx_train, ],
                       classes = "lqph")
  )
  
  oof_prob_insects_me[idx_test] <- predict(
    m_insect_fold,
    newdata = X_insects[idx_test, ],
    type = "cloglog"
  )
}

# performance metrics
eval_df_insects <- data.frame(
  observed = y_insects,
  predicted = oof_prob_insects_me
) %>% filter(!is.na(predicted))

auc_insects_me <- pROC::auc(eval_df_insects$observed, eval_df_insects$predicted)

thresh_insects_me <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:nrow(eval_df_insects), eval_df_insects),
  opt.methods="MaxSens+Spec"
)[1,"predicted"]

cm_insects_me <- PresenceAbsence::cmx(
  data.frame(id=1:nrow(eval_df_insects), eval_df_insects),
  threshold=thresh_insects_me
)

tss_insects_me <- PresenceAbsence::sensitivity(cm_insects_me) +
  PresenceAbsence::specificity(cm_insects_me) - 1

# Run MaxEnt – plants ------------------------------------------------------

# prep data
X_plants <- plant_me_df[, -1]
y_plants <- plant_me_df$presence

# initialize storage for out-of-fold predictions
oof_prob_plants_me <- rep(NA_real_, nrow(plant_me_df))

# loop over folds
for (i in seq_along(folds_plants$folds_list)) {
  idx_train <- folds_plants$folds_list[[i]][[1]]
  idx_test  <- folds_plants$folds_list[[i]][[2]]

  m_plant_fold <- maxnet::maxnet(
    p = y_plants[idx_train],
    data = X_plants[idx_train, ],
    f = maxnet::maxnet.formula(
      p = y_plants[idx_train],
      data = X_plants[idx_train, ],
      classes = "lqph"
    )
  )

  oof_prob_plants_me[idx_test] <- predict(
    m_plant_fold,
    newdata = X_plants[idx_test, ],
    type = "cloglog"
  )
}

# performance metrics
eval_df_plants <- data.frame(
  observed = y_plants,
  predicted = oof_prob_plants_me
) %>%
  dplyr::filter(!is.na(predicted))

auc_plants_me <- pROC::auc(eval_df_plants$observed, eval_df_plants$predicted)

thresh_plants_me <- PresenceAbsence::optimal.thresholds(
  data.frame(id = 1:nrow(eval_df_plants), eval_df_plants),
  opt.methods = "MaxSens+Spec"
)[1, "predicted"]

cm_plants_me <- PresenceAbsence::cmx(
  data.frame(id = 1:nrow(eval_df_plants), eval_df_plants),
  threshold = thresh_plants_me
)

tss_plants_me <- PresenceAbsence::sensitivity(cm_plants_me) +
  PresenceAbsence::specificity(cm_plants_me) - 1

# Fit final GAMs ----------------------------------------------------------

# 2-part (hurdle) GAM for insects
f_bin_insects <- as.formula(paste(
  "present ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_insects), collapse=" + ")
))
gam_bin_insects <- mgcv::gam(f_bin_insects, family=binomial(),
                             data=dat_insects, method="REML")

dat_insects_pos <- dat_insects[dat_insects$present == 1, ]
f_mu_insects <- as.formula(paste(
  "species_richness ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_insects), collapse=" + ")
))
gam_mu_insects <- mgcv::gam(f_mu_insects, family=gaussian(),
                            data=dat_insects_pos, method="REML")

# 2-part GAM for plants
f_bin_plants <- as.formula(paste(
  "present ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_plants), collapse=" + ")
))
gam_bin_plants <- mgcv::gam(f_bin_plants, family=binomial(),
                            data=dat_plants, method="REML")

dat_plants_pos <- dat_plants[dat_plants$present == 1, ]
f_mu_plants <- as.formula(paste(
  "species_richness ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_plants), collapse=" + ")
))
gam_mu_plants <- mgcv::gam(f_mu_plants, family=gaussian(),
                           data=dat_plants_pos, method="REML")

# Predict w/GAMs ----------------------------------------------------------

prob_insects <- terra::predict(pred, gam_bin_insects, type="response")
prob_plants  <- terra::predict(pred, gam_bin_plants,  type="response")
mu_insects   <- terra::predict(pred, gam_mu_insects,  type="response")
mu_plants    <- terra::predict(pred, gam_mu_plants,   type="response")

# combine for expected richness maps
mu_insects[mu_insects < 0] <- 0
expected_insects <- prob_insects * mu_insects
mu_plants[mu_plants < 0] <- 0
expected_plants <- prob_plants * mu_plants

# Save GAM predictions
writeRaster(prob_insects,  "raster/insect_prob_gam.tif", overwrite = TRUE)
writeRaster(mu_insects,    "raster/insect_rich_gam.tif", overwrite = TRUE)
writeRaster(expected_insects, "raster/insect_pred_rich_gam.tif", overwrite = TRUE)

writeRaster(prob_plants,   "raster/plant_prob_gam.tif", overwrite = TRUE)
writeRaster(mu_plants,     "raster/plant_rich_gam.tif", overwrite = TRUE)
writeRaster(expected_plants, "raster/plant_pred_rich_gam.tif", overwrite = TRUE)


# Fit final BRTs ----------------------------------------------------------

# insects
df_occ_insects <- data.frame(present = dat_insects$present,
                             dat_insects[, preds_insects, drop=FALSE])

brt_bin_insects <- gbm::gbm(
  formula = present ~ .,
  data = df_occ_insects,
  distribution = "bernoulli",
  n.trees = 5000, interaction.depth = 3,
  shrinkage = 0.01, bag.fraction = 0.6,
  cv.folds = 5, verbose = FALSE
)
best_iter_bin_insects <- gbm::gbm.perf(brt_bin_insects, method="cv",
                                       plot.it=FALSE)

dat_insects_pos <- dat_insects[dat_insects$present == 1, ]
df_mu_insects <- data.frame(species_richness = dat_insects_pos$species_richness,
                            dat_insects_pos[, preds_insects, drop=FALSE])

brt_mu_insects <- gbm::gbm(
  formula = species_richness ~ .,
  data = df_mu_insects,
  distribution = "gaussian",
  n.trees = 5000, interaction.depth = 3,
  shrinkage = 0.01, bag.fraction = 0.6,
  cv.folds = 5, verbose = FALSE
)
best_iter_mu_insects <- gbm::gbm.perf(brt_mu_insects, method="cv", 
                                      plot.it=FALSE)

# plants
df_occ_plants <- data.frame(present = dat_plants$present,
                            dat_plants[, preds_plants, drop=FALSE])

brt_bin_plants <- gbm::gbm(
  formula = present ~ .,
  data = df_occ_plants,
  distribution = "bernoulli",
  n.trees = 5000, interaction.depth = 2, 
  shrinkage = 0.01, bag.fraction = 0.8,
  cv.folds = 5, verbose = FALSE
)
best_iter_bin_plants <- gbm::gbm.perf(brt_bin_plants, method="cv", 
                                      plot.it=FALSE)

dat_plants_pos <- dat_plants[dat_plants$present == 1, ]
df_mu_plants <- data.frame(species_richness = dat_plants_pos$species_richness,
                           dat_plants_pos[, preds_plants, drop=FALSE])

brt_mu_plants <- gbm::gbm(
  formula = species_richness ~ .,
  data = df_mu_plants,
  distribution = "gaussian",
  n.trees = 5000, interaction.depth = 2,
  shrinkage = 0.01, bag.fraction = 0.8,
  cv.folds = 5, verbose = FALSE
)
best_iter_mu_plants <- gbm::gbm.perf(brt_mu_plants, method="cv", 
                                     plot.it=FALSE)

# Predict w/BRTs ----------------------------------------------------------

# insects
prob_insects_brt <- terra::predict(pred, brt_bin_insects,
  fun=function(model, data){predict.gbm(model, newdata=data,
                                        n.trees=best_iter_bin_insects,
                                        type="response")})
mu_insects_brt <- terra::predict(pred, brt_mu_insects,
  fun=function(model, data){pmax(predict.gbm(model, newdata=data,
                                             n.trees=best_iter_mu_insects,
                                             type="response"),0)})


# mask to study area
prob_insects_brt <- terra::mask(prob_insects_brt, expected_insects)
mu_insects_brt   <- terra::mask(mu_insects_brt,  expected_insects)

# combine for expected richness
expected_insects_brt <- prob_insects_brt * mu_insects_brt
plot(expected_insects_brt)

# plants
prob_plants_brt <- terra::predict(pred, brt_bin_plants,
  fun=function(model, data){predict.gbm(model, newdata=data,
                                        n.trees=best_iter_bin_plants,
                                        type="response")})
mu_plants_brt <- terra::predict(pred, brt_mu_plants,
  fun=function(model, data){pmax(predict.gbm(model, newdata=data,
                                             n.trees=best_iter_mu_plants,
                                             type="response"),0)})

prob_plants_brt <- terra::mask(prob_plants_brt, expected_plants)
mu_plants_brt   <- terra::mask(mu_plants_brt,  expected_plants)

expected_plants_brt <- prob_plants_brt * mu_plants_brt
plot(expected_plants_brt)

# Save BRT predictions
writeRaster(prob_insects_brt,  "raster/insect_prob_brt.tif", overwrite = TRUE)
writeRaster(mu_insects_brt,    "raster/insect_rich_brt.tif", overwrite = TRUE)
writeRaster(expected_insects_brt, "raster/insect_pred_rich_brt.tif", overwrite = TRUE)

writeRaster(prob_plants_brt,   "raster/plant_prob_brt.tif", overwrite = TRUE)
writeRaster(mu_plants_brt,     "raster/plant_rich_brt.tif", overwrite = TRUE)
writeRaster(expected_plants_brt, "raster/plant_pred_rich_brt.tif", overwrite = TRUE)

# Fit final RFs -----------------------------------------------------------

# insects
rf_bin_insects <- randomForest(
  x = dat_insects[, preds_insects, drop=FALSE],
  y = as.factor(dat_insects$present),
  ntree = 10000, mtry = floor(sqrt(length(preds_insects)))
)

dat_insects_pos <- dat_insects[dat_insects$present == 1, ]
rf_mu_insects <- randomForest(
  x = dat_insects_pos[, preds_insects, drop=FALSE],
  y = dat_insects_pos$species_richness,
  ntree = 10000, mtry = floor(sqrt(length(preds_insects)))
)

# plants
rf_bin_plants <- randomForest(
  x = dat_plants[, preds_plants, drop=FALSE],
  y = as.factor(dat_plants$present),
  ntree = 10000, mtry = floor(sqrt(length(preds_plants)))
)

dat_plants_pos <- dat_plants[dat_plants$present == 1, ]
rf_mu_plants <- randomForest(
  x = dat_plants_pos[, preds_plants, drop=FALSE],
  y = dat_plants_pos$species_richness,
  ntree = 10000, mtry = floor(sqrt(length(preds_plants)))
)

# Predict from RFs --------------------------------------------------------

# NB: parallelizing because it otherwise crashes

n_cores <- 20
terraOptions(
  threads = n_cores,
  tempdir = "/home/NINA.NO/jenny.hansen/Mounts/scratch/tmp_jenny" 
)

write_compact <- function(x, filename) {
  writeRaster(x, filename,
              overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=LZW", "TILED=YES"))
}

# insects
prob_insects_rf <- terra::predict(
  pred,
  rf_bin_insects,
  fun = function(model, data) {
    res <- rep(NA_real_, nrow(data))
    valid <- complete.cases(data)
    if (any(valid)) {
      probs <- predict(model,
                       newdata = as.data.frame(data[valid, , drop = FALSE]),
                       type = "prob")[, "1"]
      res[valid] <- probs
    }
    res
  },
  filename = "raster/insect_prob_rf_expanded.tif",
  overwrite = TRUE
)

mu_insects_rf <- terra::predict(
  pred,
  rf_mu_insects,
  fun = function(model, data) {
    res <- rep(NA_real_, nrow(data))
    valid <- complete.cases(data)
    if (any(valid)) {
      vals <- predict(model, newdata = as.data.frame(data[valid, , drop = FALSE]))
      res[valid] <- pmax(vals, 0)
    }
    res
  },
  filename = "raster/insect_rich_rf_expanded.tif",
  overwrite = TRUE
)

expected_insects_rf <- mask(prob_insects_rf * mu_insects_rf, expected_insects)
writeRaster(expected_insects_rf, "raster/insect_pred_rich_rf_expanded.tif", 
            overwrite = TRUE)


# plants
prob_plants_rf <- terra::predict(
  pred,
  rf_bin_plants,
  fun = function(model, data) {
    res <- rep(NA_real_, nrow(data))
    valid <- complete.cases(data)
    if (any(valid)) {
      probs <- predict(model,
                       newdata = as.data.frame(data[valid, , drop = FALSE]),
                       type = "prob")[, "1"]
      res[valid] <- probs
    }
    res
  },
  filename = "raster/plant_prob_rf.tif",
  overwrite = TRUE
)

mu_plants_rf <- terra::predict(
  pred,
  rf_mu_plants,
  fun = function(model, data) {
    res <- rep(NA_real_, nrow(data))
    valid <- complete.cases(data)
    if (any(valid)) {
      vals <- predict(model, newdata = as.data.frame(data[valid, , drop = FALSE]))
      res[valid] <- pmax(vals, 0)
    }
    res
  },
  filename = "raster/plant_rich_rf.tif",
  overwrite = TRUE
)

expected_plants_rf  <- mask(prob_plants_rf * mu_plants_rf, expected_plants)
writeRaster(expected_plants_rf,  "raster/plant_pred_rich_rf.tif", overwrite = TRUE)

# Fit final MaxEnt --------------------------------------------------------

# insects
X_insects <- insect_me_df[, -1]  
y_insects <- insect_me_df$presence  

final_maxent_insects <- maxnet(
  p = y_insects,
  data = X_insects,
  f = maxnet.formula(
    p = y_insects,
    data = X_insects,
    classes = "lqph"  
  ),
  regmult = 1.5  
)

# plants
X_plants <- plant_me_df[, -1]  
y_plants <- plant_me_df$presence  

final_maxent_plants <- maxnet(
  p = y_plants,
  data = X_plants,
  f = maxnet.formula(
    p = y_plants,
    data = X_plants,
    classes = "lqph"
  ),
  regmult = 1.5
)

# Predict from MaxEnt -----------------------------------------------------

# insects
maxnet_pred_fun <- function(model, data) {
  data_df <- as.data.frame(data)
  preds <- predict(model, newdata = data_df, type = "cloglog")
  as.numeric(preds)  # ensure vector output
}

prob_insects_me <- terra::predict(
  pred,
  final_maxent_insects,
  fun = maxnet_pred_fun,
  na.rm = TRUE
)

# mask to study area
prob_insects_me <- terra::mask(prob_insects_me, expected_insects)

# plants
prob_plants_me <- terra::predict(
  pred,
  final_maxent_plants,
  fun = maxnet_pred_fun,
  na.rm = TRUE
)

prob_plants_me <- terra::mask(prob_plants_me, expected_plants)

# Save MaxEnt predictions
writeRaster(prob_insects_me, "raster/insect_prob_maxent.tif", overwrite = TRUE)
writeRaster(prob_plants_me,  "raster/plant_prob_maxent.tif", overwrite = TRUE)


# Compute ensemble weights ------------------------------------------------

# NB: weights are based on cross-validation prediction performance, not
# final fit model

# Occurrence weights (AUC - 0.5 so no-skill = 0)
norm2 <- function(v) v / sum(v, na.rm=TRUE)

# insects 

# get conditional rmse on presence-only data
is_pos_insects <- dat_insects$present == 1

rmse_mu_insects_gam <- sqrt(mean((oof_mu_insects[is_pos_insects] -
                                    dat_insects$species_richness[is_pos_insects])^2, 
                                 na.rm=TRUE))
rmse_mu_insects_brt <- sqrt(mean((oof_mu_brt_insects[is_pos_insects] -
                                    dat_insects$species_richness[is_pos_insects])^2, 
                                 na.rm=TRUE))
rmse_mu_insects_rf  <- sqrt(mean((oof_mu_rf_insects[is_pos_insects] -
                                    dat_insects$species_richness[is_pos_insects])^2, 
                                 na.rm=TRUE))

w_occ_insects <- norm2(c(
  GAM = as.numeric(auc_occ_insects) - 0.5,
  BRT = as.numeric(auc_occ_brt_insects) - 0.5,
  RF  = as.numeric(auc_rf_insects) - 0.5,
  ME  = as.numeric(auc_insects_me) - 0.5
))

w_mu_insects <- norm2(c(
  GAM = 1/(rmse_mu_insects_gam^2),
  BRT = 1/(rmse_mu_insects_brt^2),
  RF  = 1/(rmse_mu_insects_rf^2)
))

# plants

is_pos_plants <- dat_plants$present == 1

rmse_mu_plants_gam <- sqrt(mean((oof_mu_plants[is_pos_plants] -
                                   dat_plants$species_richness[is_pos_plants])^2, 
                                na.rm=TRUE))
rmse_mu_plants_brt <- sqrt(mean((oof_mu_brt_plants[is_pos_plants] -
                                   dat_plants$species_richness[is_pos_plants])^2,
                                na.rm=TRUE))
rmse_mu_plants_rf  <- sqrt(mean((oof_mu_rf_plants[is_pos_plants] -
                                   dat_plants$species_richness[is_pos_plants])^2,
                                na.rm=TRUE))

w_occ_plants <- norm2(c(
  GAM = as.numeric(auc_occ_plants) - 0.5,
  BRT = as.numeric(auc_occ_brt_plants) - 0.5,
  RF  = as.numeric(auc_rf_plants) - 0.5,
  ME  = as.numeric(auc_plants_me) - 0.5
))

w_mu_plants <- norm2(c(
  GAM = 1/(rmse_mu_plants_gam^2),
  BRT = 1/(rmse_mu_plants_brt^2),
  RF  = 1/(rmse_mu_plants_rf^2)
))

# Ensemble model ----------------------------------------------------------

# insects

# OOF ensemble (cannot include maxent here because it has more bg pts)
ens_prob_insects <- w_occ_insects["GAM"] * oof_prob_insects +
  w_occ_insects["BRT"] * oof_prob_brt_insects +
  w_occ_insects["RF"]  * oof_prob_rf_insects

ens_mu_insects <- w_mu_insects["GAM"] * oof_mu_insects +
  w_mu_insects["BRT"] * oof_mu_brt_insects +
  w_mu_insects["RF"]  * oof_mu_rf_insects
ens_mu_insects[is.na(ens_mu_insects)] <- 0

# OOF expected richness
ensemble_pred_insects <- ens_prob_insects * ens_mu_insects

# metrics 
obs_insects <- dat_insects$species_richness
rmse_ens_insects <- sqrt(mean((ensemble_pred_insects - obs_insects)^2, na.rm=TRUE))
r2_ens_insects   <- 1 - sum((obs_insects - ensemble_pred_insects)^2, na.rm=TRUE) /
  sum((obs_insects - mean(obs_insects, na.rm=TRUE))^2, na.rm=TRUE)

auc_ens_insects <- pROC::auc(dat_insects$present, ens_prob_insects)

thresh_ens_insects <- PresenceAbsence::optimal.thresholds(
  data.frame(id=seq_len(nrow(dat_insects)),
             observed=dat_insects$present,
             predicted=ens_prob_insects),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm_ens_insects <- PresenceAbsence::cmx(
  data.frame(id=seq_len(nrow(dat_insects)),
             observed=dat_insects$present,
             predicted=ens_prob_insects),
  threshold=thresh_ens_insects)

sens_ens_insects <- PresenceAbsence::sensitivity(cm_ens_insects, st.dev=FALSE)
spec_ens_insects <- PresenceAbsence::specificity(cm_ens_insects, st.dev=FALSE)
tss_ens_insects  <- as.numeric(sens_ens_insects + spec_ens_insects - 1)

# plants

# OOF ensemble
ens_prob_plants <- w_occ_plants["GAM"] * oof_prob_plants +
  w_occ_plants["BRT"] * oof_prob_brt_plants +
  w_occ_plants["RF"]  * oof_prob_rf_plants

ens_mu_plants <- w_mu_plants["GAM"] * oof_mu_plants +
  w_mu_plants["BRT"] * oof_mu_brt_plants +
  w_mu_plants["RF"]  * oof_mu_rf_plants
ens_mu_plants[is.na(ens_mu_plants)] <- 0

ensemble_pred_plants <- ens_prob_plants * ens_mu_plants

obs_plants <- dat_plants$species_richness
rmse_ens_plants <- sqrt(mean((ensemble_pred_plants - obs_plants)^2, na.rm=TRUE))
r2_ens_plants   <- 1 - sum((obs_plants - ensemble_pred_plants)^2, na.rm=TRUE) /
  sum((obs_plants - mean(obs_plants, na.rm=TRUE))^2, na.rm=TRUE)

auc_ens_plants <- pROC::auc(dat_plants$present, ens_prob_plants)

thresh_ens_plants <- PresenceAbsence::optimal.thresholds(
  data.frame(id=seq_len(nrow(dat_plants)),
             observed=dat_plants$present,
             predicted=ens_prob_plants),
  opt.methods="MaxSens+Spec")[1,"predicted"]

cm_ens_plants <- PresenceAbsence::cmx(
  data.frame(id=seq_len(nrow(dat_plants)),
             observed=dat_plants$present,
             predicted=ens_prob_plants),
  threshold=thresh_ens_plants)

sens_ens_plants <- PresenceAbsence::sensitivity(cm_ens_plants, st.dev=FALSE)
spec_ens_plants <- PresenceAbsence::specificity(cm_ens_plants, st.dev=FALSE)
tss_ens_plants  <- as.numeric(sens_ens_plants + spec_ens_plants - 1)


# Predict from ensemble model ---------------------------------------------

# insects
prob_insects_ens <- prob_insects * w_occ_insects["GAM"] +
  prob_insects_brt  * w_occ_insects["BRT"] +
  prob_insects_rf   * w_occ_insects["RF"]  +
  prob_insects_me   * w_occ_insects["ME"]   


mu_insects_ens <- w_mu_insects["GAM"] * mu_insects +
  w_mu_insects["BRT"] * mu_insects_brt +
  w_mu_insects["RF"]  * mu_insects_rf
mu_insects_ens <- terra::clamp(mu_insects_ens, lower = 0, upper = Inf)

expected_insects_ens <- prob_insects_ens * mu_insects_ens
expected_insects_ens <- terra::mask(expected_insects_ens, expected_insects)
plot(expected_insects_ens)
writeRaster(expected_insects_ens, "raster/ensemble_prediction_insects.tif",
            overwrite = TRUE)

# plants
prob_plants_ens <- prob_plants * w_occ_plants["GAM"] +
  prob_plants_brt  * w_occ_plants["BRT"] +
  prob_plants_rf   * w_occ_plants["RF"]  +
  prob_plants_me   * w_occ_plants["ME"]   

mu_plants_ens <- w_mu_plants["GAM"] * mu_plants +
  w_mu_plants["BRT"] * mu_plants_brt +
  w_mu_plants["RF"]  * mu_plants_rf
mu_plants_ens <- terra::clamp(mu_plants_ens, lower = 0, upper = Inf)

expected_plants_ens <- prob_plants_ens * mu_plants_ens
expected_plants_ens <- terra::mask(expected_plants_ens, expected_plants)
plot(expected_plants_ens)
writeRaster(expected_plants_ens, "raster/ensemble_prediction_plants.tif",
            overwrite = TRUE)

# Comparison table --------------------------------------------------------

# Convert all TSS objects to plain numeric values
tss_occ_insects     <- as.numeric(tss_occ_insects[1])
tss_occ_brt_insects <- as.numeric(tss_occ_brt_insects[1])
tss_rf_insects      <- as.numeric(tss_rf_insects[1])
tss_ens_insects     <- as.numeric(tss_ens_insects[1])
tss_insects_me      <- as.numeric(tss_insects_me[1]) 

tss_occ_plants      <- as.numeric(tss_occ_plants[1])
tss_occ_brt_plants  <- as.numeric(tss_occ_brt_plants[1])
tss_rf_plants       <- as.numeric(tss_rf_plants[1])
tss_ens_plants      <- as.numeric(tss_ens_plants[1])
tss_plants_me       <- as.numeric(tss_plants_me[1])   


# summary table 
summary_df <- tibble::tribble(
  ~Taxon,   ~Model,     ~RMSE,                ~R2,                 ~AUC,                          ~TSS,
  # insects
  "Insects","GAM",      rmse_exp_insects,     r2_exp_insects,      as.numeric(auc_occ_insects),   tss_occ_insects,
  "Insects","BRT",      rmse_brt_insects,     r2_brt_insects,      as.numeric(auc_occ_brt_insects), tss_occ_brt_insects,
  "Insects","RF",       rmse_rf_insects,      r2_rf_insects,       as.numeric(auc_rf_insects),    tss_rf_insects,
  "Insects","MaxEnt",   NA,                   NA,                  as.numeric(auc_insects_me),    tss_insects_me,  # <- added
  "Insects","Ensemble", rmse_ens_insects,     r2_ens_insects,      as.numeric(auc_ens_insects),   tss_ens_insects,
  
  # plants
  "Plants", "GAM",      rmse_exp_plants,      r2_exp_plants,       as.numeric(auc_occ_plants),    tss_occ_plants,
  "Plants", "BRT",      rmse_brt_plants,      r2_brt_plants,       as.numeric(auc_occ_brt_plants), tss_occ_brt_plants,
  "Plants", "RF",       rmse_rf_plants,       r2_rf_plants,        as.numeric(auc_rf_plants),     tss_rf_plants,
  "Plants", "MaxEnt",   NA,                   NA,                  as.numeric(auc_plants_me),     tss_plants_me,   # <- added
  "Plants", "Ensemble", rmse_ens_plants,      r2_ens_plants,       as.numeric(auc_ens_plants),    tss_ens_plants
) %>%
  dplyr::mutate(across(where(is.numeric), \(x) round(x, 3)))


summary_df
saveRDS(summary_df, "models/model_comparison_table.rds")

# Comparison plot ---------------------------------------------------------

# pivot data
summary_long <- summary_df %>%
  tidyr::pivot_longer(cols = c(RMSE, R2, AUC, TSS),
                      names_to = "Metric", values_to = "Value")

# ensure model order 
summary_long$Model <- factor(summary_long$Model, 
                             levels = c("GAM", "BRT", "RF", "MaxEnt", "Ensemble"))

# compute relative performance within each Taxon × Metric group
summary_long <- summary_long %>%
  group_by(Taxon, Metric) %>%
  mutate(
    # invert RMSE so lower = better
    Score = ifelse(Metric == "RMSE", -Value, Value),
    RelScore = (Score - min(Score, na.rm = TRUE)) /
      (max(Score, na.rm = TRUE) - min(Score, na.rm = TRUE))
  ) %>%
  ungroup()

# RandomForest by far the best performing models for insects
# Boosted Regression Trees performed best for plants, but
# the ensemble model was close

# plot showing performance
ggplot(summary_long, aes(x = Metric, y = Value, 
                         fill = RelScore, group = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Taxon, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(title = "Model performance comparison",
       y = "Metric value", x = "Performance metric", fill = "Relative\nperformance") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Model),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3, color = "black", angle = 45, hjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) 


# Export models -----------------------------------------------------------

# save individual models 

# insects
saveRDS(gam_bin_insects, "models/gam_bin_insects.rds")
saveRDS(gam_mu_insects,  "models/gam_mu_insects.rds")
saveRDS(brt_bin_insects, "models/brt_bin_insects.rds")
saveRDS(brt_mu_insects,  "models/brt_mu_insects.rds")
saveRDS(rf_bin_insects,  "models/rf_bin_insects.rds")
saveRDS(rf_mu_insects,   "models/rf_mu_insects.rds")
saveRDS(ens_mu_insects,  "models/ens_mu_insects.rds")
saveRDS(final_maxent_insects, "models/maxent_insects.rds")  

# plants
saveRDS(gam_bin_plants, "models/gam_bin_plants.rds")
saveRDS(gam_mu_plants,  "models/gam_mu_plants.rds")
saveRDS(brt_bin_plants, "models/brt_bin_plants.rds")
saveRDS(brt_mu_plants,  "models/brt_mu_plants.rds")
saveRDS(rf_bin_plants,  "models/rf_bin_plants.rds")
saveRDS(rf_mu_plants,   "models/rf_mu_plants.rds")
saveRDS(ens_mu_plants,  "models/ens_mu_plants.rds")
saveRDS(final_maxent_plants, "models/maxent_plants.rds")    

# all models in one list 

all_models <- list(
  insects = list(
    GAM = list(bin = gam_bin_insects, mu = gam_mu_insects),
    BRT = list(bin = brt_bin_insects, mu = brt_mu_insects),
    RF  = list(bin = rf_bin_insects,  mu = rf_mu_insects),
    ME  = final_maxent_insects, 
    ENS = ens_mu_insects
  ),
  plants = list(
    GAM = list(bin = gam_bin_plants, mu = gam_mu_plants),
    BRT = list(bin = brt_bin_plants, mu = brt_mu_plants),
    RF  = list(bin = rf_bin_plants,  mu = rf_mu_plants),
    ME  = final_maxent_plants,    
    ENS = ens_mu_plants
  )
)

saveRDS(all_models, "models/all_models.rds")

