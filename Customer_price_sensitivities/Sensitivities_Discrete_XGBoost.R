##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Affiliation:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         18/02/2022

##### Clear workspace #####
cat("\014")
rm(list = ls())
while (dev.cur()>1) dev.off()

#### Libraries #####
Packages <- c('dplyr', 'fastDummies', 'Matrix', 'lme4', 'gbm', 'xgboost', 'glmnet')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load output from matching procedure
# Load results from Matching_Discrete.R

##### Matching #####
### Matching settings
Match <- matrix(list(), N, T_max)
# Number of closest matches to find
Prop_Rank <- 10
# Determine realized offers
Realized_i <- sapply(1:N, function(i) which(Z[i, ] == 1))
# Find observations that did receive offer t
N_t <- sapply(1:T_max, function(t) which(Z[, t] == 1))

### Match indices
# Determine the matches for every policy i
for (i in 1:N) {
  # Match for every counterfactual rate change category
  for (t in (1:T_max)[-Realized_i[i]]) {
    # Compute the absolute propensity score distance
    Distance_prob <- abs(Prop_Score[i, t] - Prop_Score[N_t[[t]], t])
    # Compute the ranked distances
    Rank_Distance <- rank(Distance_prob, ties.method = 'first')
    # Find the matches that are closest in distance
    Match[i, t][[1]] <- AMI$Churn[N_t[[t]][which(Rank_Distance <= Prop_Rank)]]
  }
}

### Multiple imputation
# Settings
M <- 10
Augm_Response <- array(NA, dim = c(N, T_max, M))
# Imputation
for (m in 1:M) {
  Augm_Response[, , m][which(Z == 1)] <- AMI$Churn
}
# Find the corresponding counterfactual offers
Counterfactual_i <- t(sapply(1:N, function(i) (1:T_max)[-Realized_i[i]]))

### Counterfactual sampling
# Ensure reproducibility by setting the seed for the random number generator
seed <- 2019
set.seed(seed)
for (i in 1:N) {
  # Sample from matched responses to infer counterfactual response
  Augm_Response[i, Counterfactual_i[i, ], ] <- t(sapply(Counterfactual_i[i, ], function(t) 
    rbinom(n = M, size = 1, prob = mean(Match[i, t][[1]]))))
}

### Augmented dataset
AMI_Augm <- as.data.frame(lapply(AMI, rep, each = T_max, byrow = TRUE))
AMI_Augm$Bucket <- as.vector(rep(Q_buckets, times = N))
AMI_Augm$Offer <- as.vector(rep(c(1:T_max), times = N))
AMI_Augm$Offer <- relevel(factor(AMI_Augm$Offer), ref = names(which.max(table(AMI$Offer))))

### Average churn
Avg_Churn <- as.vector(sapply(1:N, function(i) sapply(1:T_max, function(t) sum(Augm_Response[i, t, ]) / M)))
AMI_Augm$Churn <- as.integer(ifelse(Avg_Churn >= 0.5, 1, 0))

##### Causal inference #####
### Input variables
Var_XGB <- as.matrix(dummy_cols(AMI_Augm[, c('Churn', 'Competitiveness', 'Offer',
                                             'Policy_Type', 'Risk_Group', 'Undershooting_1', 
                                             'Undershooting_2', 'Premium_New_Base')], 
                                remove_most_frequent_dummy = FALSE)[, -c(1, 3, 4, 5)])
Var_XGB <- Matrix(Var_XGB, sparse = TRUE)
Labels_XGB <- as.vector(AMI_Augm$Churn)
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)

### Settings for Extreme Gradient Boosting (XGBoost), with hyperparameters optimized 
### for the continuous approach in Sensitivities_Continuous.R
# General settings
Churn_Trees <- 10000
Churn_Stop <- 250
seed <- 2019
# Optimal grid
Churn_Grid_Opt <- expand.grid(
  eta = 0.03,
  max_depth = 10,
  min_child_weight = 0,
  subsample = 0.9,
  colsample_bytree = 0.6,
  gamma = 0,
  lambda = 0,
  alpha = 0,
  logloss = NA,
  Trees = NA
)

### Given the optimized hyperparameters, determine optimal number of decision trees
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Define the parameters of the XGBoost model
par_XGB <- list(objective = 'binary:logistic', eval_metric = 'logloss', 
                booster = 'gbtree', seed = seed,
                eta = Churn_Grid_Opt$eta, 
                max_depth = Churn_Grid_Opt$max_depth, 
                min_child_weight = Churn_Grid_Opt$min_child_weight,
                subsample = Churn_Grid_Opt$subsample, 
                colsample_bytree = Churn_Grid_Opt$colsample_bytree,
                gamma = Churn_Grid_Opt$gamma, lambda = Churn_Grid_Opt$lambda, 
                alpha = Churn_Grid_Opt$alpha)
# Perform cross-validation on the XGBoost model
Churn_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                    nfold = 10, nrounds = Churn_Trees, early_stopping_rounds = Churn_Stop, 
                    verbose = TRUE)
# Store the minimal ASAM obtained and the corresponding number of decision trees used
Churn_Grid_Opt$logloss <- min(Churn_XGB$evaluation_log$test_logloss_mean)
Churn_Grid_Opt$Trees <- which.min(Churn_XGB$evaluation_log$test_logloss_mean)

### Form the cross-validated XGBoost model
# Optimal number of decision trees
Churn_Best <- Churn_Grid_Opt$Trees
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Define the parameters of the XGBoost model
par_XGB <- list(objective = 'binary:logistic', eval_metric = 'logloss',
                booster = 'gbtree', seed = seed,
                eta = Churn_Grid_Opt$eta, 
                max_depth = Churn_Grid_Opt$max_depth, 
                min_child_weight = Churn_Grid_Opt$min_child_weight,
                subsample = Churn_Grid_Opt$subsample, 
                colsample_bytree = Churn_Grid_Opt$colsample_bytree,
                gamma = Churn_Grid_Opt$gamma, lambda = Churn_Grid_Opt$lambda, 
                alpha = Churn_Grid_Opt$alpha)
# Estimate the cross-validated XGBoost model
Churn_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                     nrounds = Churn_Best,
                     verbose = TRUE)

### Churn predictions
Churn_Pred <- predict(Churn_XGB, newdata = Var_XGB, type = 'response')
Churn_Pred_matrix <- t(sapply(1:N, function(i) Churn_Pred[((i - 1) * T_max + 1):(i * T_max)]))

##### Average expected churn #####
### Churn prediction function
Churn_Pred <- function(Churn_Model, Offers, Comp, Data_Row) {
  # Competitiveness adjustment due to interdependencies
  if (!is.na(Comp)) {
    Data_Row$Competitiveness <- Comp
    Data_Row$Premium_New_Base <- Data_Row$Risk_Base * Data_Row$D_1 / (Comp + 1)
    Data_Row$Undershooting_1 <- pmax(Data_Row$D_1 - Data_Row$Premium_New_Base, 0, na.rm = TRUE)
    Data_Row$Undershooting_2 <-pmax(Data_Row$D_2 - Data_Row$Premium_New_Base, 0, na.rm = TRUE)
  }
  
  # Offers to consider
  if (nrow(Data_Row) == 1) {
    Data_Row <- as.data.frame(lapply(Data_Row, rep, each = length(Offers), byrow = TRUE))
  }
  Data_Row$Offer[] <- Offers
  
  # Estimated churn rate
  Var_Row_XGB <- as.matrix(dummy_cols(Data_Row[, c('Churn', 'Competitiveness', 'Offer',
                                                   'Policy_Type', 'Risk_Group', 'Undershooting_1', 
                                                   'Undershooting_2', 'Premium_New_Base')], 
                                      remove_most_frequent_dummy = FALSE)[, -c(1, 3, 4, 5)])
  Var_Row_XGB <- Matrix(Var_Row_XGB, sparse = TRUE)
  Churn <- predict(Churn_Model, newdata = as.matrix(Var_Row_XGB), type = 'response')
  
  # Return the output
  return(Churn)
}

### Price sensitivity surface with competitiveness
# Preliminary variables
T_Offers <- sort(unique(AMI$Rate_Change), decreasing = FALSE)
T_grid <- seq(from = ceiling(min(T_Offers) * 100) / 100, to = floor(max(T_Offers) * 100) / 100, 
              by = 0.01)
T_grid <- c(min(AMI$Price_Change), T_grid, max(AMI$Price_Change))
T_Comp <- sort(unique(AMI$Competitiveness), decreasing = FALSE)
C_grid <- seq(from = ceiling(min(T_Comp) * 100) / 100, to = floor(max(T_Comp) * 100) / 100, 
              by = 0.02)
C_grid <- c(min(AMI$Competitiveness), C_grid, max(AMI$Competitiveness))
C_max <- length(C_grid)
Comp_GRM <- matrix(NA, T_max, C_max)
# Aggregatation
for (t in 1:T_max) {
  # Average expected churn
  Comp_GRM[t, ] <- 1 / N * sapply(1:C_max, function(c)
    sum(Churn_Pred(Churn_Model = Churn_XGB,
                   Offers = as.character(t), Comp = C_grid[c],
                   Data_Row = AMI)))
}

### Price sensitivity surface with level of risk
# Preliminary variables
R_max <- nlevels(AMI$Risk_Group)
R_levels <- levels(AMI$Risk_Group)
Risk_Offer_GRM <- matrix(NA, R_max, T_max)
AMI_r <- AMI
# Aggregation
for (r in 1:R_max) {
  # Level of risk under consideration
  AMI_r$Risk_Group[] <- R_levels[r]
  # Average expected churn
  Risk_Offer_GRM[r, ] <- 1 / N * sapply(1:T_max, function(t)
    sum(Churn_Pred(Churn_Model = Churn_XGB, Offers = as.character(t), 
                   Comp = NA, Data_Row = AMI_r)))
}

### Average surface for level of risk and competitiveness
# Preliminary variables
R_max <- nlevels(AMI$Risk_Group)
R_levels <- levels(AMI$Risk_Group)
Risk_Comp_GRM <- matrix(NA, R_max, C_max)
AMI_r <- AMI
# Aggregation
for (r in 1:R_max) {
  # Level of risk under consideration
  AMI_r$Risk_Group[] <- R_levels[r]
  # Average expected churn
  Risk_Comp_GRM[r, ] <- 1 / N * sapply(1:C_max, function(c)
    sum(Churn_Pred(Churn_Model = Churn_XGB, Offers = AMI_r$Offer,
                   Comp = C_grid[c], Data_Row = AMI_r)))
}

### Store the final results
