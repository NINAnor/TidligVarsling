# script to select predictor variables for models
# Jenny Hansen
# 27 September 2025

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(tidyverse)
library(randomForestSRC)

# Import data -------------------------------------------------------------

# data prepared in script 13
insect_df <- read.csv("data/insect_model_df.csv") 
plant_df  <- read.csv("data/plant_model_df.csv") 


# Identify correlated vars ------------------------------------------------

# function to determine variable correlation
select_uncorrelated_joint <- function(df, importance_df, cutoff = 0.7) {
  
  cor_mat <- cor(df, use = "pairwise.complete.obs")
  
  ranked_vars <- importance_df %>%
    arrange(desc(Importance_joint)) %>%
    pull(Variable)
  
  excluded <- c()
  summary_list <- list()
  cluster_id <- 1
  
  for (var in ranked_vars) {
    if (var %in% excluded) next
    
    correlated <- names(which(abs(cor_mat[var, ]) >= cutoff))
    correlated <- setdiff(correlated, var)
    
    summary_list[[cluster_id]] <- tibble(
      Cluster = cluster_id,
      Kept    = var,
      Dropped = ifelse(length(correlated) == 0, NA,
                       paste(correlated, collapse = ", "))
    )
    
    excluded <- c(excluded, correlated)
    cluster_id <- cluster_id + 1
  }
  
  bind_rows(summary_list)
}


# Joint RF models ---------------------------------------------------------

# function to fit joint randomforest models
joint_rf_importance <- function(plant_df, insect_df,
                                       ntree = 5000,
                                       replicates = 10,
                                       nsplit = 5,
                                       combine_fun = c("mean", "max")) {
  combine_fun <- match.arg(combine_fun)
  
  # join datasets
  # need to 'pad out' plants with NA to match nrow for insects
  plant_vec <- rep(NA_real_, nrow(insect_df))
  plant_vec[seq_len(nrow(plant_df))] <- plant_df$species_richness
  
  joint_df <- insect_df %>%
    select(-species_richness) %>%
    mutate(
      insect = insect_df$species_richness,
      plant  = plant_vec
    )
  
  
  # tune mtry and nodesize
  message("Tuning mtry and nodesize...")
  tune_res <- tune(
    Multivar(plant, insect) ~ .,
    data = joint_df,
    ntreeTry = 10000
  )
  
  opt_mtry     <- tune_res$optimal[["mtry"]]
  opt_nodesize <- tune_res$optimal[["nodesize"]]
  message(glue::glue("Optimal mtry = {opt_mtry}, nodesize = {opt_nodesize}"))
  
  # replicate forests
  set.seed(123)
  imp_list <- replicate(replicates, {
    fit <- rfsrc(
      Multivar(plant, insect) ~ .,
      data = joint_df,
      importance = "permute",
      ntree = ntree,
      mtry = opt_mtry,
      nodesize = opt_nodesize,
      nsplit = nsplit,
      samptype = "swor" # try without replacement for stability
    )
    
    tibble(
      Variable = names(fit$regrOutput$plant$importance),
      Plant    = fit$regrOutput$plant$importance,
      Insect   = fit$regrOutput$insect$importance
    )
  }, simplify = FALSE)
  
  # aggregate importance across replicates
  imp_stable <- bind_rows(imp_list, .id = "replicate") %>%
    group_by(Variable) %>%
    summarise(
      Importance_plant  = mean(Plant),
      Importance_insect = mean(Insect),
      .groups = "drop"
    )
  
  
  # normalize per group
  imp_stable <- imp_stable %>%
    mutate(
      Plant_norm  = Importance_plant  / max(Importance_plant, na.rm = TRUE),
      Insect_norm = Importance_insect / max(Importance_insect, na.rm = TRUE),
      Importance_joint = (Plant_norm + Insect_norm) / 2, 
      Importance_joint_max = pmax(Plant_norm, Insect_norm)
    )
  
  return(imp_stable %>% arrange(desc(Importance_joint)))
}


# run the models
imp_joint <- joint_rf_importance(
  plant_df, insect_df,
  ntree = 10000,
  replicates = 10
) %>%
  select(Variable,
         Importance_plant,
         Importance_insect,
         Importance_joint,
         Importance_joint_max) # Optimal mtry = 9, nodesize = 1

imp_joint

# importance values for plants have absolutely tanked when including
# only potentially new invasives

# model selection
shared_selection_summary <- select_uncorrelated_joint(
  df = plant_df %>% select(-species_richness),
  importance_df = imp_joint,
  cutoff = 0.7) 

# Selection table ---------------------------------------------------------

kept_vars <- shared_selection_summary %>% pull(Kept)

final_table <- imp_joint %>%
  filter(Variable %in% kept_vars) %>%
  arrange(desc(Importance_joint))

final_table_ranked <- final_table %>%
  mutate(Rank = row_number())
final_table_ranked

# pivot table for ggplot
plot_data <- final_table_ranked %>%
  pivot_longer(cols = starts_with("Importance"),
    names_to = "Group", values_to = "Importance") %>%
  mutate(Group = recode(Group,
                        Importance_plant  = "Plant",
                        Importance_insect = "Insect",
                        Importance_joint  = "Joint"))


# Plot --------------------------------------------------------------------

# unique ranks and labels from the ranked table
rank_breaks <- final_table_ranked$Rank
rank_labels <- final_table_ranked$Variable

# plot
ggplot(plot_data, aes(x = Rank, y = Importance, color = Group)) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_x_continuous(
    breaks = rank_breaks,
    labels = rank_labels
  ) +
  labs(
    title = "Variable importance drop-off",
    x = "Predictor (ranked by joint importance)",
    y = "Importance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8)
  )


# top ten vars, jointly ranked
final_table_ranked %>% print(n = 10)

# all vars, ranked
final_table_ranked %>% print(n = 18)


# Thoughts on selection ---------------------------------------------------

# original models had 14 variables- I personally like to follow the
# 1 variable per n = 10, which would be 7
# there are stronger signals for plants than insects, so I suggest using
# the ranking for plants for variable selection
# if we follow the variable importance curve for plants, it drops at 
# distance to garden center
# if we include all variables to that point, it will be 11

# NB: this has changed dramatically after the plant filtering
# suggest keeping distance to avfall, even though that means upping the
# variables to 12
