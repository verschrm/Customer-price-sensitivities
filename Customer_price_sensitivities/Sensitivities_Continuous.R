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
Packages <- c('dplyr', 'fastDummies', 'Matrix', 'lme4', 'gbm', 'xgboost')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load output from matching procedure
# Load results from Matching_Continuous.R

##### Causal inference #####
### Input variables
Var_XGB <- as.matrix(dummy_cols(AMI[, c('Churn', 'Competitiveness', 'Offer', 'GPS',
                                        'Policy_Type', 'Risk_Group', 'Undershooting_1',
                                        'Undershooting_2', 'Premium_New_Base')], 
                                remove_most_frequent_dummy = FALSE)[, -c(1, 5, 6)])
Var_XGB <- Matrix(Var_XGB, sparse = TRUE)
Labels_XGB <- as.vector(AMI$Churn)
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)

### Settings for Extreme Gradient Boosting (XGBoost)
# General settings
Churn_Trees <- 10000
Churn_Stop <- 250
seed <- 2019
# Grid search settings
Churn_eta <- c(seq(from = 0.01, to = 0.05, by = 0.01), seq(from = 0.1, to = 0.25, by = 0.05), 0.5)
Churn_Depth <- c(0, 1, seq(from = 2, to = 10, by = 2), 25, 50)
Churn_Min <- c(0, seq(from = 1, to = 5, by = 1), 10, 25, 50)
Churn_row_sample <- c(seq(from = 0.1, to = 1, by = 0.1))
Churn_col_sample <- c(seq(from = 0.1, to = 1, by = 0.1))
Churn_gamma <- c(0, 0.1, 1, 10, 100)
Churn_lambda <- c(0, 0.1, 1, 10, 100)
Churn_alpha <- c(0, 0.1, 1, 10, 100)
# Optimal grid
Churn_Grid_Opt <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = NA,
  colsample_bytree = NA,
  gamma = NA,
  lambda = NA,
  alpha = NA,
  logloss = NA,
  Trees = NA
)
# Grid 1 - Learning rate
Churn_Grid_1 <- expand.grid(
  eta = Churn_eta,
  max_depth = 6,
  min_child_weight = 1,
  subsample = 1,
  colsample_bytree = 1,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  logloss = NA,
  Trees = NA
)
# Grid 2 - Tree-specific parameters
Churn_Grid_2 <- expand.grid(
  eta = NA,
  max_depth = Churn_Depth,
  min_child_weight = Churn_Min,
  subsample = 1,
  colsample_bytree = 1,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  logloss = NA,
  Trees = NA
)
# Grid 3 - Stochastic features
Churn_Grid_3 <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = Churn_row_sample,
  colsample_bytree = Churn_col_sample,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  logloss = NA,
  Trees = NA
)
# Grid 4 - Regularization penalties
Churn_Grid_4 <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = NA,
  colsample_bytree = NA,
  gamma = Churn_gamma,
  lambda = Churn_lambda,
  alpha = Churn_alpha,
  logloss = NA,
  Trees = NA
)

### Grid search 1
for (i in 1:nrow(Churn_Grid_1)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(objective = 'binary:logistic', eval_metric = 'logloss',
                  booster = 'gbtree', seed = seed, eta = Churn_Grid_1$eta[i], 
                  max_depth = Churn_Grid_1$max_depth[i], 
                  min_child_weight = Churn_Grid_1$min_child_weight[i],
                  subsample = Churn_Grid_1$subsample[i], 
                  colsample_bytree = Churn_Grid_1$colsample_bytree[i],
                  gamma = Churn_Grid_1$gamma[i], lambda = Churn_Grid_1$lambda[i], 
                  alpha = Churn_Grid_1$alpha[i])
  # Perform cross-validation on the XGBoost model
  Churn_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                      nfold = 10, nrounds = Churn_Trees, 
                      early_stopping_rounds = Churn_Stop, verbose = TRUE)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  Churn_Grid_1$logloss[i] <- min(Churn_XGB$evaluation_log$test_logloss_mean)
  Churn_Grid_1$Trees[i] <- which.min(Churn_XGB$evaluation_log$test_logloss_mean)
}
# Store the optimal parameters from the grid
Churn_Grid_Opt$eta <- Churn_Grid_1$eta[which.min(Churn_Grid_1$logloss)]
# Adopt these parameters in the next grid searches
Churn_Grid_2$eta <- Churn_Grid_Opt$eta
Churn_Grid_3$eta <- Churn_Grid_Opt$eta
Churn_Grid_4$eta <- Churn_Grid_Opt$eta
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_Churn_2 <- which(apply(sapply(1:(ncol(Churn_Grid_Opt) - 2), function(j) 
  Churn_Grid_2[, j] == Churn_Grid_1[which.min(Churn_Grid_1$logloss), j]), 1, sum) == (ncol(Churn_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
Churn_Grid_2$logloss[Excl_Churn_2] <- Churn_Grid_1$logloss[which.min(Churn_Grid_1$logloss)]
Churn_Grid_2$Trees[Excl_Churn_2] <- Churn_Grid_1$Trees[which.min(Churn_Grid_1$logloss)]

### Grid search 2
for (i in c(1:nrow(Churn_Grid_2))[-Excl_Churn_2]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(objective = 'binary:logistic', eval_metric = 'logloss',
                  booster = 'gbtree', seed = seed, eta = Churn_Grid_2$eta[i], 
                  max_depth = Churn_Grid_2$max_depth[i], 
                  min_child_weight = Churn_Grid_2$min_child_weight[i],
                  subsample = Churn_Grid_2$subsample[i], 
                  colsample_bytree = Churn_Grid_2$colsample_bytree[i],
                  gamma = Churn_Grid_2$gamma[i], lambda = Churn_Grid_2$lambda[i], 
                  alpha = Churn_Grid_2$alpha[i])
  # Perform cross-validation on the XGBoost model
  Churn_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                    nfold = 10, nrounds = Churn_Trees, 
                    early_stopping_rounds = Churn_Stop, verbose = TRUE)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  Churn_Grid_2$logloss[i] <- min(Churn_XGB$evaluation_log$test_logloss_mean)
  Churn_Grid_2$Trees[i] <- which.min(Churn_XGB$evaluation_log$test_logloss_mean)
}
# Store the optimal parameters from the grid
Churn_Grid_Opt$max_depth <- Churn_Grid_2$max_depth[which.min(Churn_Grid_2$logloss)]
Churn_Grid_Opt$min_child_weight <- Churn_Grid_2$min_child_weight[which.min(Churn_Grid_2$logloss)]
# Adopt these parameters in the next grid searches
Churn_Grid_3$max_depth <- Churn_Grid_Opt$max_depth
Churn_Grid_3$min_child_weight <- Churn_Grid_Opt$min_child_weight
Churn_Grid_4$max_depth <- Churn_Grid_Opt$max_depth
Churn_Grid_4$min_child_weight <- Churn_Grid_Opt$min_child_weight
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_Churn_3 <- which(apply(sapply(1:(ncol(Churn_Grid_Opt) - 2), function(j) 
  Churn_Grid_3[, j] == Churn_Grid_2[which.min(Churn_Grid_2$logloss), j]), 1, sum) == (ncol(Churn_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
Churn_Grid_3$logloss[Excl_Churn_3] <- Churn_Grid_2$logloss[which.min(Churn_Grid_2$logloss)]
Churn_Grid_3$Trees[Excl_Churn_3] <- Churn_Grid_2$Trees[which.min(Churn_Grid_2$logloss)]

### Grid search 3
for (i in c(1:nrow(Churn_Grid_3))[-Excl_Churn_3]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(objective = 'binary:logistic', eval_metric = 'logloss',
                  booster = 'gbtree', seed = seed, eta = Churn_Grid_3$eta[i], 
                  max_depth = Churn_Grid_3$max_depth[i], 
                  min_child_weight = Churn_Grid_3$min_child_weight[i],
                  subsample = Churn_Grid_3$subsample[i], 
                  colsample_bytree = Churn_Grid_3$colsample_bytree[i],
                  gamma = Churn_Grid_3$gamma[i], lambda = Churn_Grid_3$lambda[i], 
                  alpha = Churn_Grid_3$alpha[i])
  # Perform cross-validation on the XGBoost model
  Churn_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                    nfold = 10, nrounds = Churn_Trees, 
                    early_stopping_rounds = Churn_Stop, verbose = TRUE)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  Churn_Grid_3$logloss[i] <- min(Churn_XGB$evaluation_log$test_logloss_mean)
  Churn_Grid_3$Trees[i] <- which.min(Churn_XGB$evaluation_log$test_logloss_mean)
}
# Store the optimal parameters from the grid
Churn_Grid_Opt$subsample <- Churn_Grid_3$subsample[which.min(Churn_Grid_3$logloss)]
Churn_Grid_Opt$colsample_bytree <- Churn_Grid_3$colsample_bytree[which.min(Churn_Grid_3$logloss)]
# Adopt these parameters in the next grid search
Churn_Grid_4$subsample <- Churn_Grid_Opt$subsample
Churn_Grid_4$colsample_bytree <- Churn_Grid_Opt$colsample_bytree
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_Churn_4 <- which(apply(sapply(1:(ncol(Churn_Grid_Opt) - 2), function(j) 
  Churn_Grid_4[, j] == Churn_Grid_3[which.min(Churn_Grid_3$logloss), j]), 1, sum) == (ncol(Churn_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
Churn_Grid_4$logloss[Excl_Churn_4] <- Churn_Grid_3$logloss[which.min(Churn_Grid_3$logloss)]
Churn_Grid_4$Trees[Excl_Churn_4] <- Churn_Grid_3$Trees[which.min(Churn_Grid_3$logloss)]

### Grid search 4
for (i in c(1:nrow(Churn_Grid_4))[-Excl_Churn_4]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(objective = 'binary:logistic', eval_metric = 'logloss', 
                  booster = 'gbtree', seed = seed, eta = Churn_Grid_4$eta[i], 
                  max_depth = Churn_Grid_4$max_depth[i], 
                  min_child_weight = Churn_Grid_4$min_child_weight[i],
                  subsample = Churn_Grid_4$subsample[i], 
                  colsample_bytree = Churn_Grid_4$colsample_bytree[i],
                  gamma = Churn_Grid_4$gamma[i], lambda = Churn_Grid_4$lambda[i], 
                  alpha = Churn_Grid_4$alpha[i])
  # Perform cross-validation on the XGBoost model
  Churn_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                      nfold = 10, nrounds = Churn_Trees, 
                      early_stopping_rounds = Churn_Stop, verbose = TRUE)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  Churn_Grid_4$logloss[i] <- min(Churn_XGB$evaluation_log$test_logloss_mean)
  Churn_Grid_4$Trees[i] <- which.min(Churn_XGB$evaluation_log$test_logloss_mean)
}
# Store the optimal parameters from the grid
Churn_Grid_Opt$gamma <- Churn_Grid_4$gamma[which.min(Churn_Grid_4$logloss)]
Churn_Grid_Opt$lambda <- Churn_Grid_4$lambda[which.min(Churn_Grid_4$logloss)]
Churn_Grid_Opt$alpha <- Churn_Grid_4$alpha[which.min(Churn_Grid_4$logloss)]
Churn_Grid_Opt$logloss <- Churn_Grid_4$logloss[which.min(Churn_Grid_4$logloss)]
Churn_Grid_Opt$Trees <- Churn_Grid_4$Trees[which.min(Churn_Grid_4$logloss)]

### Form the cross-validated XGBoost model
# Optimal number of decision trees
Churn_Best <- Churn_Grid_Opt$Trees
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Define the parameters of the XGBoost model
par_XGB <- list(objective = 'binary:logistic', eval_metric = 'logloss',
                booster = 'gbtree', seed = seed, eta = Churn_Grid_Opt$eta, 
                max_depth = Churn_Grid_Opt$max_depth, 
                min_child_weight = Churn_Grid_Opt$min_child_weight,
                subsample = Churn_Grid_Opt$subsample, 
                colsample_bytree = Churn_Grid_Opt$colsample_bytree,
                gamma = Churn_Grid_Opt$gamma, lambda = Churn_Grid_Opt$lambda, 
                alpha = Churn_Grid_Opt$alpha)
# Estimate the cross-validated XGBoost model
Churn_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                     nrounds = Churn_Best, verbose = TRUE)

### Churn predictions
# For observed treatments only, as the counterfactual treatments are a continuum of treatments
Churn_Pred_vector <- predict(Churn_XGB, newdata = as.matrix(Var_XGB), type = 'response')

##### Average expected churn #####
### Churn prediction function
Churn_Pred <- function(Churn_Model, Offers, Comp, GPS_Model, sd_GPS, sd_GPS_i, Data_Row) {
  # Competitiveness adjustment due to interdependencies
  if (!is.na(Comp)) {
    Data_Row$Competitiveness <- Comp
    Data_Row$Premium_New_Base <- Data_Row$Risk_Base * Data_Row$D_1 / (Comp + 1)
    Data_Row$Undershooting_1 <- pmax(Data_Row$D_1 - Data_Row$Premium_New_Base, 0, na.rm = TRUE)
    Data_Row$Undershooting_2 <- pmax(Data_Row$D_2 - Data_Row$Premium_New_Base, 0, na.rm = TRUE)
  }
  
  # Generalized propensity score input
  Var_Row_XGB <- as.matrix(dummy_cols(Data_Row[, c('Offer', 'Policy_Type', 'Risk_Group', 
                                                   'Undershooting_1', 'Undershooting_2', 
                                                   'Premium_New_Base')], 
                                  remove_most_frequent_dummy = FALSE)[, -c(1, 2, 3)])
  Var_Row_XGB <- Matrix(Var_Row_XGB, sparse = TRUE)
  Mean <- predict(GPS_Model, newdata = as.matrix(Var_Row_XGB), type = 'response')
  Normalizer <- apply(sapply(1:(length(Q_bins) - 1), function(t) 
    pnorm((Q_bins[t + 1] - Mean) / sd_GPS[t], mean = 0, sd = 1) - 
      pnorm((Q_bins[t] - Mean) / sd_GPS[t], mean = 0, sd = 1)), 1, sum)
  
  # Offers to consider
  if (nrow(Data_Row) == 1) {
    Data_Row <- as.data.frame(lapply(Data_Row, rep, each = length(Offers), byrow = TRUE))
  }
  Data_Row$Offer <- Offers
  
  # Estimated GPS
  Data_Row$GPS <- dnorm((Offers - Mean) / sd_GPS_i, mean = 0, sd = 1) / Normalizer / sd_GPS_i
  
  # Estimated churn rate
  Var_Row_XGB <- as.matrix(dummy_cols(Data_Row[, c('Competitiveness', 'Offer', 'GPS',
                                                   'Policy_Type', 'Risk_Group', 'Undershooting_1',
                                                   'Undershooting_2', 'Premium_New_Base')], 
                                      remove_most_frequent_dummy = FALSE)[, -c(4, 5)])
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
T_sd <- sd_GPS_XGB[as.numeric(cut(T_grid, breaks = Q_bins, include.lowest = TRUE))]
T_Comp <- sort(unique(AMI$Competitiveness), decreasing = FALSE)
C_grid <- seq(from = ceiling(min(T_Comp) * 100) / 100, to = floor(max(T_Comp) * 100) / 100, 
              by = 0.02)
C_grid <- c(min(AMI$Competitiveness), C_grid, max(AMI$Competitiveness))
C_max <- length(C_grid)
Comp_DRF <- matrix(NA, length(T_grid), C_max)
# Aggregatation
for (t in 1:length(T_grid)) {
  # Average expected churn
  Comp_DRF[t, ] <- 1 / N * sapply(1:C_max, function(c) 
    sum(Churn_Pred(Churn_Model = Churn_XGB, Offers = T_grid[t], Comp = C_grid[c], 
                   GPS_Model = GPS_XGB, sd_GPS = sd_GPS_XGB, sd_t = T_sd[t], 
                   Data_Row = AMI)))
}

### Price sensitivity surface with level of risk
# Preliminary variables
R_max <- nlevels(AMI$Risk_Group)
R_levels <- levels(AMI$Risk_Group)
Risk_Offer_DRF <- matrix(NA, R_max, length(T_grid))
AMI_r <- AMI
# Aggregation
for (r in 1:R_max) {
  # Level of risk under consideration
  AMI_r$Risk_Group[] <- R_levels[r]
  # Average expected churn
  Risk_Offer_DRF[r, ] <- 1 / N * sapply(1:T_max, function(t)
    sum(Churn_Pred(Churn_Model = Churn_XGB, Offers = T_grid[t], Comp = NA,
                   GPS_Model = GPS_XGB, sd_GPS = sd_GPS_XGB, sd_t = T_sd[t], 
                   Data_Row = AMI_r)))
}

### Average surface for level of risk and competitiveness
# Preliminary variables
R_max <- nlevels(AMI$Risk_Group)
R_levels <- levels(AMI$Risk_Group)
Risk_Comp_DRF <- matrix(NA, R_max, C_max)
AMI_r <- AMI
# Aggregation
for (r in 1:R_max) {
  # Level of risk under consideration
  AMI_r$Risk_Group[] <- R_levels[r]
  # Average expected churn
  Risk_Comp_DRF[r, ] <- 1 / N * sapply(1:C_max, function(c)
    sum(Churn_Pred(Churn_Model = Churn_XGB, Offers = AMI_r$Offer, Comp = C_grid[c],
                   GPS_Model = GPS_XGB, sd_GPS = sd_GPS_XGB, sd_t = sd_GPS_t, 
                   Data_Row = AMI_r)))
}

### Store the final results
