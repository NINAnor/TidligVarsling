# script to fit GAM and boosted regression tree models for plants and insects
# test predictive ability and generate weights for ensemble models
# Jenny Hansen
# 29 September 2025

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
library(ggplot2)
library(viridis)


# Import data -------------------------------------------------------------

# variables were selected in script 14
# NB: these changed after filtering plants and re-doing the joint RF
keep_vars <- c("species_richness",
               "neighbor_prct_impervious", "percentage_buildings", 
               "distance_to_lumberyard", "ndvi_stdDev",
               "distance_to_industrial_area", "distance_to_private_road", 
               "ndwi_median", "ndvi_summer", "neighbor_prct_open",
               "percentage_agriculture", "coldest_winter_temperature",
               "distance_to_avfall")

insect_sf <- st_read("vector/insect_model_sf.geojson")  %>% 
  select(all_of(keep_vars)) %>% 
  mutate(group = "insect",
         present = as.integer(species_richness > 0))
plant_sf  <- st_read("vector/plant_model_sf.geojson") %>% 
  select(all_of(keep_vars))  %>% 
  mutate(group = "plant",
         present = as.integer(species_richness > 0))

# predictors from script 10
pred <- rast("raster/complete_prediction_stack.tif")
pred <- pred[[keep_vars[2:13]]]

rm(keep_vars)


# Set up spatial folds ----------------------------------------------------

set.seed(42)

# determine spatial autocorrelation range (using raster stack)
cv_range <- cv_spatial_autocor(pred, num_sample = 50000, plot = TRUE)
cv_range$range # 9963.496

# using point data (insects)
cv_range_points <- cv_spatial_autocor(x = insect_sf, column = "present", 
                                      plot = TRUE)
cv_range_points$range # 3600.758

# using point data (plants)
cv_range_points <- cv_spatial_autocor(x = plant_sf, column = "present", 
                                      plot = TRUE)
cv_range_points$range # 6538.501

# changing block approach to account for smaller plant dataset
block_i = 3600
block_p <- 6500

folds_insects <- cv_spatial(x = insect_sf, "present", k = 5, size = block_i,
                            selection = "random", iteration = 200)
folds_plants <- cv_spatial(x = plant_sf, "present", k = 5, size = block_p,
                           selection = "random", iteration = 200)


# GAM hurdle- insects -----------------------------------------------------

# insects

# set up data
preds_insects <- names(insect_sf %>% 
                 select(neighbor_prct_impervious:distance_to_avfall) %>% 
                 st_drop_geometry())
dat_insects <- insect_sf %>% st_drop_geometry()
dat_insects$present <- as.integer(dat_insects$species_richness > 0)

# out-of-fold storage
oof_prob_insects <- rep(NA_real_, nrow(dat_insects))  # probability of >0
oof_mu_insects   <- rep(NA_real_, nrow(dat_insects))  # conditional mean richness


# 2-part hurdle model 
for (i in seq_along(folds_insects$folds_list)) {
  idx_train <- folds_insects$folds_list[[i]][[1]]  # training indices
  idx_test  <- folds_insects$folds_list[[i]][[2]]  # test indices
  
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

rmse_exp_insects
r2_exp_insects

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

auc_occ_insects
tss_occ_insects

# GAM hurdle- plants ------------------------------------------------------

# set up data
preds_plants <- names(plant_sf %>% 
                         select(neighbor_prct_impervious:distance_to_avfall) %>% 
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
  
  # occurrence model (binomial gam)
  f_bin <- as.formula(paste("present ~",
                            paste(sprintf("s(%s,k=4,bs='cs')", preds_plants), 
                                  collapse=" + ")))
  gam_bin <- mgcv::gam(f_bin, family=binomial(), data=train, method="REML")
  
  oof_prob_plants[idx_test] <- predict(gam_bin, newdata=test, type="response")
  
  # conditional richness model (Gaussian gam)
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

# metrics on the conditional part
rmse_exp_plants <- sqrt(mean((expected_richness_plants - obs_plants)^2, 
                              na.rm=TRUE))
r2_exp_plants   <- 1 - sum((obs_plants - expected_richness_plants)^2, 
                            na.rm=TRUE) /
  sum((obs_plants - mean(obs_plants, na.rm=TRUE))^2, na.rm=TRUE)

rmse_exp_plants
r2_exp_plants # much lower R2 after filtering

# metrics on occurrence part
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

auc_occ_plants
tss_occ_plants


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

# metrics on conditional part
expected_richness_brt_insects <- oof_prob_brt_insects * oof_mu_brt_insects
rmse_brt_insects <- sqrt(mean((expected_richness_brt_insects - 
                                 dat_insects$species_richness)^2, 
                      na.rm=TRUE))
r2_brt_insects <- 1 - sum((dat_insects$species_richness - 
                             expected_richness_brt_insects)^2, 
                    na.rm=TRUE) /
  sum((dat_insects$species_richness - mean(dat_insects$species_richness))^2, 
      na.rm=TRUE)

rmse_brt_insects
r2_brt_insects

# metrics on occurrence part
y_occ_brt_insects <- dat_insects$present
prob_occ_brt_insects <- oof_prob_brt_insects

# AUC
auc_occ_brt_insects <- pROC::auc(y_occ_brt_insects, 
                                 prob_occ_brt_insects)

# optimal threshold (maximising Sens+Spec)
thresh_brt_insects <- PresenceAbsence::optimal.thresholds(
  data.frame(id=1:length(y_occ_brt_insects),
             observed=y_occ_brt_insects,
             predicted=prob_occ_brt_insects),
  opt.methods="MaxSens+Spec")[1,"predicted"]

# confusion matrix at opt. threshold
cm_brt_insects <- PresenceAbsence::cmx(
  data.frame(id=1:length(y_occ_brt_insects),
             observed=y_occ_brt_insects,
             predicted=prob_occ_brt_insects),
  threshold=thresh_brt_insects)

# TSS = sensitivity + specificity – 1
tss_occ_brt_insects <- PresenceAbsence::sensitivity(cm_brt_insects) +
  PresenceAbsence::specificity(cm_brt_insects) - 1

auc_occ_brt_insects
tss_occ_brt_insects


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
  
  # occurrence model (bernoulli) 
  df_occ <- data.frame(present = train$present, train[, preds_plants, 
                                                      drop=FALSE])
  
  # adjusted due to smaller sample size
  brt_bin <- gbm::gbm(
    formula = present ~ .,
    data = df_occ,
    distribution = "bernoulli",
    n.trees = 5000, interaction.depth = 2,
    shrinkage = 0.01, bag.fraction = 0.8,
    n.minobsinnode = 5,
    cv.folds = 5, verbose = FALSE
  )
  best_iter_bin <- gbm::gbm.perf(brt_bin, method="cv", plot.it=FALSE)
  
  oof_prob_brt_plants[idx_test] <- predict(
    brt_bin,
    newdata = test[, preds_plants, drop=FALSE],
    n.trees = best_iter_bin,
    type = "response"
  )
  
  # conditional richness model (Gaussian)
  train_pos <- train[train$present == 1, ]
  if (nrow(train_pos) >= 10) {
    df_mu <- data.frame(species_richness = train_pos$species_richness,
                        train_pos[, preds_plants, drop=FALSE])
    
    brt_mu <- gbm::gbm(
      formula = species_richness ~ .,
      data = df_mu,
      distribution = "gaussian",
      n.trees = 5000, interaction.depth = 2,
      shrinkage = 0.01, bag.fraction = 0.8,
      cv.folds = 5, verbose = FALSE
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

# metrics for conditional part
expected_richness_brt_plants <- oof_prob_brt_plants * oof_mu_brt_plants
rmse_brt_plants <- sqrt(mean((expected_richness_brt_plants - 
                                dat_plants$species_richness)^2, 
                      na.rm=TRUE))
r2_brt_plants   <- 1 - sum((dat_plants$species_richness - 
                              expected_richness_brt_plants)^2, 
                    na.rm=TRUE) /
  sum((dat_plants$species_richness - mean(dat_plants$species_richness))^2, 
      na.rm=TRUE)

rmse_brt_plants
r2_brt_plants

# metrics for occurrence part
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

auc_occ_brt_plants
tss_occ_brt_plants



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

# metrics (insects)
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


rmse_rf_insects
r2_rf_insects
auc_rf_insects
tss_rf_insects


# RF hurdle-plants --------------------------------------------------------

oof_prob_rf_plants <- rep(NA_real_, nrow(dat_plants))
oof_mu_rf_plants   <- rep(NA_real_, nrow(dat_plants))

for (i in seq_along(folds_plants$folds_list)) {
  idx_train <- folds_plants$folds_list[[i]][[1]]
  idx_test  <- folds_plants$folds_list[[i]][[2]]
  
  train <- dat_plants[idx_train, ]
  test  <- dat_plants[idx_test, ]
  
  # occurrence model
  rf_bin <- randomForest(
    x = train[, preds_plants, drop=FALSE],
    y = as.factor(train$present),
    ntree = 1000, mtry = floor(sqrt(length(preds_plants)))
  )
  
  oof_prob_rf_plants[idx_test] <- predict(
    rf_bin, newdata = test[, preds_plants, drop=FALSE],
    type = "prob")[, "1"]
  
  # conditional richness model
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

# metrics (plants)
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

rmse_rf_plants
r2_rf_plants
auc_rf_plants
tss_rf_plants


# Fit final GAMs ----------------------------------------------------------

# 2-part (hurdle) GAM for insects

# occurrence (binomial)
f_bin_insects <- as.formula(paste(
  "present ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_insects), collapse=" + ")
))
gam_bin_insects <- mgcv::gam(f_bin_insects, family=binomial(),
                             data=dat_insects, method="REML")

# conditional richness (Gaussian)
dat_insects_pos <- dat_insects[dat_insects$present == 1, ]
f_mu_insects <- as.formula(paste(
  "species_richness ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_insects), collapse=" + ")
))
gam_mu_insects <- mgcv::gam(f_mu_insects, family=gaussian(),
                            data=dat_insects_pos, method="REML")


summary(gam_mu_insects)
plot(gam_mu_insects)


# 2-part GAM for plants

# occurrence
f_bin_plants <- as.formula(paste(
  "present ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_plants), collapse=" + ")
))
gam_bin_plants <- mgcv::gam(f_bin_plants, family=binomial(),
                            data=dat_plants, method="REML")

# conditional richness
dat_plants_pos <- dat_plants[dat_plants$present == 1, ]
f_mu_plants <- as.formula(paste(
  "species_richness ~",
  paste(sprintf("s(%s,k=4,bs='cs')", preds_plants), collapse=" + ")
))
gam_mu_plants <- mgcv::gam(f_mu_plants, family=gaussian(),
                           data=dat_plants_pos, method="REML")

summary(gam_mu_plants)
plot(gam_mu_plants) # oh boy- this is a weird one


# Predict w/GAMs ----------------------------------------------------------

# predict occurrence probability 
prob_insects <- terra::predict(pred, gam_bin_insects, type="response")
prob_plants  <- terra::predict(pred, gam_bin_plants,  type="response")

# predict conditional richness 
mu_insects <- terra::predict(pred, gam_mu_insects, type="response")
mu_plants  <- terra::predict(pred, gam_mu_plants,  type="response")

# combine for expected richness maps
mu_insects[mu_insects < 0] <- 0
expected_insects <- prob_insects * mu_insects
mu_plants[mu_plants < 0] <- 0
expected_plants <- prob_plants * mu_plants


plot(expected_insects)
plot(expected_plants)


# Fit final BRTs ----------------------------------------------------------

# insects

# prob of presence
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

# conditional richness
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

summary(brt_mu_insects)

# plants

# prob of presence
df_occ_plants <- data.frame(present = dat_plants$present,
                             dat_plants[, preds_plants, drop=FALSE])

brt_bin_plants <- gbm::gbm(
  formula = present ~ .,
  data = df_occ_plants,
  distribution = "bernoulli",
  n.trees = 5000, interaction.depth = 2, # adjusted for smaller n
  shrinkage = 0.01, bag.fraction = 0.8,
  cv.folds = 5, verbose = FALSE
)
best_iter_bin_plants <- gbm::gbm.perf(brt_bin_plants, method="cv", 
                                      plot.it=FALSE)

# conditional richness
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

summary(brt_mu_plants)


# Predict w/BRTs ----------------------------------------------------------

# insects

# predict occurrence probability 
prob_insects_brt <- terra::predict(
  pred, 
  brt_bin_insects,
  fun = function(model, data) {
    predict.gbm(model, newdata = data,
                n.trees = best_iter_bin_insects,
                type = "response")
  }
)

# predict conditional richness (clamped to 0)
mu_insects_brt <- terra::predict(
  pred,
  brt_mu_insects,
  fun = function(model, data) {
    pred_vals <- predict.gbm(model, newdata = data,
                             n.trees = best_iter_mu_insects,
                             type = "response")
    pmax(pred_vals, 0)  # clamp negatives
  }
)


# mask to study area
prob_insects_brt <- terra::mask(prob_insects_brt, expected_insects)
mu_insects_brt   <- terra::mask(mu_insects_brt,  expected_insects)

# combine for expected richness
expected_insects_brt <- prob_insects_brt * mu_insects_brt
plot(expected_insects_brt)

# plants

# predict occurrence probability 
prob_plants_brt <- terra::predict(
  pred, 
  brt_bin_plants,
  fun = function(model, data) {
    predict.gbm(model, newdata = data,
                n.trees = best_iter_bin_plants,
                type = "response")
  }
)

# predict conditional richness (clamped to 0)
mu_plants_brt <- terra::predict(
  pred,
  brt_mu_plants,
  fun = function(model, data) {
    pred_vals <- predict.gbm(model, newdata = data,
                             n.trees = best_iter_mu_plants,
                             type = "response")
    pmax(pred_vals, 0)  
  }
)


# mask to study area
prob_plants_brt <- terra::mask(prob_plants_brt, expected_plants)
mu_plants_brt   <- terra::mask(mu_plants_brt,  expected_plants)

# combine for expected richness
expected_plants_brt <- prob_plants_brt * mu_plants_brt
plot(expected_plants_brt)



# Fit final RFs -----------------------------------------------------------

# insects
# occurrence (binomial)
rf_bin_insects <- randomForest(
  x = dat_insects[, preds_insects, drop=FALSE],
  y = as.factor(dat_insects$present),
  ntree = 1000, mtry = floor(sqrt(length(preds_insects)))
)

# conditional richness (regression, presence only)
dat_insects_pos <- dat_insects[dat_insects$present == 1, ]
rf_mu_insects <- randomForest(
  x = dat_insects_pos[, preds_insects, drop=FALSE],
  y = dat_insects_pos$species_richness,
  ntree = 1000, mtry = floor(sqrt(length(preds_insects)))
)

# plants
# occurrence (binomial)
rf_bin_plants <- randomForest(
  x = dat_plants[, preds_plants, drop=FALSE],
  y = as.factor(dat_plants$present),
  ntree = 1000, mtry = floor(sqrt(length(preds_plants)))
)

# conditional richness (regression, presence only)
dat_plants_pos <- dat_plants[dat_plants$present == 1, ]
rf_mu_plants <- randomForest(
  x = dat_plants_pos[, preds_plants, drop=FALSE],
  y = dat_plants_pos$species_richness,
  ntree = 1000, mtry = floor(sqrt(length(preds_plants)))
)


# Predict from RFs --------------------------------------------------------

# insects
# predict occurrence probability
prob_insects_rf <- terra::predict(
  pred,
  rf_bin_insects,
  fun = function(model, data) {
    apply(data, 1, function(row) {
      if (any(is.na(row))) {
        return(NA_real_)
      } else {
        predict(model, newdata = as.data.frame(t(row)), type = "prob")[, "1"]
      }
    })
  }
)


# predict conditional richness (clamp to 0)
mu_insects_rf <- terra::predict(
  pred,
  rf_mu_insects,
  fun = function(model, data) {
    apply(data, 1, function(row) {
      if (any(is.na(row))) {
        return(NA_real_)
      } else {
        val <- predict(model, newdata = as.data.frame(t(row)))
        max(val, 0)  # clamp negatives
      }
    })
  }
)

# mask to study area
prob_insects_rf <- terra::mask(prob_insects_rf, expected_insects)
mu_insects_rf   <- terra::mask(mu_insects_rf, expected_insects)

# combine for expected richness
expected_insects_rf <- prob_insects_rf * mu_insects_rf
plot(expected_insects_rf)


# plants

# predict occurrence probability
prob_plants_rf <- terra::predict(
  pred,
  rf_bin_plants,
  fun = function(model, data) {
    apply(data, 1, function(row) {
      if (any(is.na(row))) {
        return(NA_real_)
      } else {
        predict(model, newdata = as.data.frame(t(row)), type = "prob")[, "1"]
      }
    })
  }
)


# predict conditional richness (clamp to 0)
mu_plants_rf <- terra::predict(
  pred,
  rf_mu_plants,
  fun = function(model, data) {
    apply(data, 1, function(row) {
      if (any(is.na(row))) {
        return(NA_real_)
      } else {
        val <- predict(model, newdata = as.data.frame(t(row)))
        max(val, 0)  # clamp negatives
      }
    })
  }
)

# mask to study area
prob_plants_rf <- terra::mask(prob_plants_rf, expected_plants)
mu_plants_rf   <- terra::mask(mu_plants_rf, expected_plants)

# combine for expected richness
expected_plants_rf <- prob_plants_rf * mu_plants_rf
plot(expected_plants_rf)

# Compute ensemble weights ------------------------------------------------

# weights are based on cross-validation prediction performance, not
# final fit model

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

# Occurrence weights (AUC - 0.5 so no-skill = 0)
norm2 <- function(v) v / sum(v, na.rm=TRUE)

w_occ_insects <- norm2(c(
  GAM = as.numeric(auc_occ_insects) - 0.5,
  BRT = as.numeric(auc_occ_brt_insects) - 0.5,
  RF  = as.numeric(auc_rf_insects) - 0.5
))

w_mu_insects <- norm2(c(
  GAM = 1/(rmse_mu_insects_gam^2),
  BRT = 1/(rmse_mu_insects_brt^2),
  RF  = 1/(rmse_mu_insects_rf^2)
))

w_occ_insects; w_mu_insects


# plants

# conditional rmse
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
  RF  = as.numeric(auc_rf_plants) - 0.5
))

w_mu_plants <- norm2(c(
  GAM = 1/(rmse_mu_plants_gam^2),
  BRT = 1/(rmse_mu_plants_brt^2),
  RF  = 1/(rmse_mu_plants_rf^2)
))

w_occ_plants; w_mu_plants



# Ensemble model ----------------------------------------------------------

# insects

# OOF ensemble
ens_prob_insects <- w_occ_insects["GAM"] * oof_prob_insects +
  w_occ_insects["BRT"] * oof_prob_brt_insects +
  w_occ_insects["RF"]  * oof_prob_rf_insects

ens_mu_insects <- w_mu_insects["GAM"] * oof_mu_insects +
  w_mu_insects["BRT"] * oof_mu_brt_insects +
  w_mu_insects["RF"]  * oof_mu_rf_insects
ens_mu_insects[is.na(ens_mu_insects)] <- 0


# OOF expected richness
ensemble_pred_insects <- ens_prob_insects * ens_mu_insects

# metrics (continuous)
obs_insects <- dat_insects$species_richness
rmse_ens_insects <- sqrt(mean((ensemble_pred_insects - obs_insects)^2, na.rm=TRUE))
r2_ens_insects   <- 1 - sum((obs_insects - ensemble_pred_insects)^2, na.rm=TRUE) /
  sum((obs_insects - mean(obs_insects, na.rm=TRUE))^2, na.rm=TRUE)

# metrics (occurrence)
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


rmse_ens_insects
r2_ens_insects
auc_ens_insects
tss_ens_insects


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

rmse_ens_plants
r2_ens_plants
auc_ens_plants
tss_ens_plants


# Predict from ensemble model ---------------------------------------------

# insects
prob_insects_ens <- w_occ_insects["GAM"] * prob_insects +
  w_occ_insects["BRT"] * prob_insects_brt +
  w_occ_insects["RF"]  * prob_insects_rf

mu_insects_ens <- w_mu_insects["GAM"] * mu_insects +
  w_mu_insects["BRT"] * mu_insects_brt +
  w_mu_insects["RF"]  * mu_insects_rf
mu_insects_ens <- terra::clamp(mu_insects_ens, lower = 0, upper = Inf)

expected_insects_ens <- prob_insects_ens * mu_insects_ens
expected_insects_ens <- terra::mask(expected_insects_ens, expected_insects)
plot(expected_insects_ens)

# plants
prob_plants_ens <- w_occ_plants["GAM"] * prob_plants +
  w_occ_plants["BRT"] * prob_plants_brt +
  w_occ_plants["RF"]  * prob_plants_rf

mu_plants_ens <- w_mu_plants["GAM"] * mu_plants +
  w_mu_plants["BRT"] * mu_plants_brt +
  w_mu_plants["RF"]  * mu_plants_rf
mu_plants_ens <- terra::clamp(mu_plants_ens, lower = 0, upper = Inf)

expected_plants_ens <- prob_plants_ens * mu_plants_ens
expected_plants_ens <- terra::mask(expected_plants_ens, expected_plants)
plot(expected_plants_ens)


# Comparison table --------------------------------------------------------

# Convert all TSS objects to plain numeric values
tss_occ_insects <- as.numeric(tss_occ_insects[1])
tss_occ_brt_insects <- as.numeric(tss_occ_brt_insects[1])
tss_rf_insects <- as.numeric(tss_rf_insects[1])
tss_ens_insects <- as.numeric(tss_ens_insects[1])

tss_occ_plants <- as.numeric(tss_occ_plants[1])
tss_occ_brt_plants <- as.numeric(tss_occ_brt_plants[1])
tss_rf_plants <- as.numeric(tss_rf_plants[1])
tss_ens_plants <- as.numeric(tss_ens_plants[1])


# summary table
summary_df <- tibble::tribble(
  ~Taxon,   ~Model,     ~RMSE,                ~R2,                 ~AUC,                        ~TSS,
  "Insects","GAM",      rmse_exp_insects,     r2_exp_insects,      as.numeric(auc_occ_insects), tss_occ_insects,
  "Insects","BRT",      rmse_brt_insects,     r2_brt_insects,      as.numeric(auc_occ_brt_insects), tss_occ_brt_insects,
  "Insects","RF",       rmse_rf_insects,      r2_rf_insects,       as.numeric(auc_rf_insects),  tss_rf_insects,
  "Insects","Ensemble", rmse_ens_insects,     r2_ens_insects,      as.numeric(auc_ens_insects), tss_ens_insects,
  "Plants", "GAM",      rmse_exp_plants,      r2_exp_plants,       as.numeric(auc_occ_plants),  tss_occ_plants,
  "Plants", "BRT",      rmse_brt_plants,      r2_brt_plants,       as.numeric(auc_occ_brt_plants), tss_occ_brt_plants,
  "Plants", "RF",       rmse_rf_plants,       r2_rf_plants,        as.numeric(auc_rf_plants),   tss_rf_plants,
  "Plants", "Ensemble", rmse_ens_plants,      r2_ens_plants,       as.numeric(auc_ens_plants),  tss_ens_plants
) %>%
  dplyr::mutate(across(where(is.numeric), \(x) round(x, 3)))


summary_df


# Comparison plot ---------------------------------------------------------

# pivot data
summary_long <- summary_df %>%
  tidyr::pivot_longer(cols = c(RMSE, R2, AUC, TSS),
               names_to = "Metric", values_to = "Value")

# ensure Model order 
summary_long$Model <- factor(summary_long$Model, 
                             levels = c("GAM","BRT","RF","Ensemble"))

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

# plot showing performance with viridis colots
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

# individual models
saveRDS(gam_bin_insects, "models/gam_bin_insects.rds")
saveRDS(gam_mu_insects,  "models/gam_mu_insects.rds")
saveRDS(brt_bin_insects, "models/brt_bin_insects.rds")
saveRDS(brt_mu_insects,  "models/brt_mu_insects.rds")
saveRDS(rf_bin_insects,  "models/rf_bin_insects.rds")
saveRDS(rf_mu_insects,   "models/rf_mu_insects.rds")
saveRDS(ens_mu_insects,  "models/ens_mu_insects.rds")

saveRDS(gam_bin_plants, "models/gam_bin_plants.rds")
saveRDS(gam_mu_plants,  "models/gam_mu_plants.rds")
saveRDS(brt_bin_plants, "models/brt_bin_plants.rds")
saveRDS(brt_mu_plants,  "models/brt_mu_plants.rds")
saveRDS(rf_bin_plants,  "models/rf_bin_plants.rds")
saveRDS(rf_mu_plants,   "models/rf_mu_plants.rds")
saveRDS(ens_mu_plants,  "models/ens_mu_plants.rds")


# all in a list
all_models <- list(
  insects = list(
    GAM = list(bin = gam_bin_insects, mu = gam_mu_insects),
    BRT = list(bin = brt_bin_insects, mu = brt_mu_insects),
    RF  = list(bin = rf_bin_insects,  mu = rf_mu_insects),
    ENS = ens_mu_insects
  ),
  plants = list(
    GAM = list(bin = gam_bin_plants, mu = gam_mu_plants),
    BRT = list(bin = brt_bin_plants, mu = brt_mu_plants),
    RF  = list(bin = rf_bin_plants,  mu = rf_mu_plants),
    ENS = ens_mu_plants
  )
)

saveRDS(all_models, "models/all_models.rds")

