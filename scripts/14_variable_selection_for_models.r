# script to select predictor variables for models
# Jenny Hansen
# 27 September 2025
# updated 08 October 2025 to include new survey data
# NB: after discussion, we decided to abandon variable
# selection and use all predictors in the models

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(tidyverse)
library(randomForest)

# Import data -------------------------------------------------------------

# data prepared in script 13
insect_df <- read.csv("data/insect_model_df.csv") 
plant_df  <- read.csv("data/plant_model_df.csv") %>% 
  select(-matches("^species_richness\\.1$"), -source)


# Identify correlated vars ------------------------------------------------

# function to determine variable correlation
select_uncorrelated_joint <- function(df, importance_df, cutoff = 0.7) {
  
  cor_mat <- cor(df, use = "pairwise.complete.obs")
  
  ranked_vars <- importance_df %>%
    filter(Variable %in% colnames(cor_mat)) %>%       
    arrange(desc(Importance_joint)) %>%
    pull(Variable)
  
  excluded <- c()
  summary_list <- list()
  cluster_id <- 1
  
  for (var in ranked_vars) {
    if (!(var %in% colnames(cor_mat))) next
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


# RF models ---------------------------------------------------------------

# function to fit individual randomforest models
rf_importance_single <- function(df, response, ntree = 500, replicates = 3) {
  vars <- setdiff(names(df), response)
  
  # drop rows with NA in predictors or response
  df_clean <- df %>%
    select(all_of(c(response, vars))) %>%
    na.omit()
  
  message(glue::glue(" - {nrow(df) - nrow(df_clean)} rows removed due to missing values ({response})"))
  
  model_list <- list()
  imp_list <- list()
  
  for (i in seq_len(replicates)) {
    set.seed(i)
    fit <- randomForest(
      reformulate(vars, response),
      data = df_clean,
      importance = TRUE,
      ntree = ntree
    )
    model_list[[i]] <- fit
    
    imp_list[[i]] <- tibble(
      Variable = vars,
      Importance = importance(fit)[, 1]
    )
  }
  
  # average importance across replicates
  imp_stable <- bind_rows(imp_list, .id = "replicate") %>%
    group_by(Variable) %>%
    summarise(Importance = mean(Importance, na.rm = TRUE), .groups = "drop")
  
  return(list(importance = imp_stable, models = model_list))
}


# function to combine importance from models
rf_importance_joint <- function(plant_df, insect_df, ntree = 500, 
                                replicates = 3) {
  
  message("Running plant random forest...")
  plant_out <- rf_importance_single(plant_df, response = "species_richness",
                                    ntree = ntree, replicates = replicates)
  imp_plant <- plant_out$importance %>% rename(Importance_plant = Importance)
  
  message("Running insect random forest...")
  insect_out <- rf_importance_single(insect_df, response = "species_richness",
                                     ntree = ntree, replicates = replicates)
  imp_insect <- insect_out$importance %>% rename(Importance_insect = Importance)
  
  # Combine & normalize importance
  imp_joint <- full_join(imp_plant, imp_insect, by = "Variable") %>%
    mutate(
      Plant_norm  = Importance_plant  / max(Importance_plant, na.rm = TRUE),
      Insect_norm = Importance_insect / max(Importance_insect, na.rm = TRUE),
      Importance_joint = (Plant_norm + Insect_norm) / 2,
      Importance_joint_max = pmax(Plant_norm, Insect_norm, na.rm = TRUE)
    ) %>%
    arrange(desc(Importance_joint))
  
  # return both importance & models
  return(list(
    importance = imp_joint,
    plant_models = plant_out$models,
    insect_models = insect_out$models,
    imp_plant = imp_plant,
    imp_insect = imp_insect
  ))
}


# run the models
rf_joint_out <- rf_importance_joint(
  plant_df, insect_df,
  ntree = 1000,
  replicates = 3
)

# extract importance
imp_plant  <- rf_joint_out$imp_plant
imp_insect <- rf_joint_out$imp_insect

imp_joint <- full_join(imp_plant, imp_insect, by = "Variable") %>%
  mutate(
    Plant_z  = (Importance_plant  - mean(Importance_plant,  na.rm = TRUE)) / sd(Importance_plant,  na.rm = TRUE),
    Insect_z = (Importance_insect - mean(Importance_insect, na.rm = TRUE)) / sd(Importance_insect, na.rm = TRUE),
    Importance_joint = (Plant_z + Insect_z) / 2,
    Importance_joint_max = pmax(Plant_z, Insect_z, na.rm = TRUE)
  ) %>%
  arrange(desc(Importance_joint))


# individual models
plant_models  <- rf_joint_out$plant_models
insect_models <- rf_joint_out$insect_models

imp_joint


# model selection
shared_selection_summary <- select_uncorrelated_joint(
  df = plant_df %>%
    sf::st_drop_geometry() %>%
    select(where(is.numeric), -species_richness) %>%
    select(all_of(intersect(names(.), imp_joint$Variable))),
  importance_df = imp_joint,
  cutoff = 0.7
)

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
final_table_ranked %>% print(n = 20)
