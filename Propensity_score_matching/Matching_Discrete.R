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
### Load confidential insurance data on automobile insurance
#   - AMI containing information on the policy renewals.

##### Outlier removal #####
### Disregard outliers in extreme rate changes and competitiveness levels
summary(AMI$Rate_Change)
Q1_Rate_Change <- quantile(AMI$Rate_Change, 0.25)
Q3_Rate_Change <- quantile(AMI$Rate_Change, 0.75)
IQR_Rate_Change <- Q3_Rate_Change - Q1_Rate_Change
summary(AMI$Competitiveness)
Q1_Comp <- quantile(AMI$Competitiveness, 0.25)
Q3_Comp <- quantile(AMI$Competitiveness, 0.75)
IQR_Comp <- Q3_Comp - Q1_Comp
AMI <- AMI[which((AMI$Rate_Change >= (Q1_Rate_Change - 1.5 * IQR_Rate_Change)) &
                   (AMI$Rate_Change <= (Q3_Rate_Change + 1.5 * IQR_Rate_Change))), ]
AMI <- AMI[which((AMI$Competitiveness >= (Q1_Comp - 1.5 * IQR_Comp)) &
                   (AMI$Competitiveness <= (Q3_Comp + 1.5 * IQR_Comp))), ]

### Resulting reduced sample size
N <- nrow(AMI)

### Buckets for discrete treatment categories
# 5 buckets
Q_grid <- 1 / 5
Q_bins <- quantile(AMI$Rate_Change, seq(0, 1, Q_grid))
T_max <- length(Q_bins) - 1

### Missing counterfactual indicator
Z <- Matrix(0, N, T_max, sparse = TRUE)
AMI$Bucket <- cut(AMI$Rate_Change, breaks = Q_bins, include.lowest = TRUE)
Q_buckets <- levels(AMI$Bucket)
AMI$Bucket <- relevel(AMI$Bucket, ref = names(which.max(table(AMI$Bucket))))
for (t in 1:T_max) {
  Z[, t] <- 1 * (Q_buckets[t] == AMI$Bucket)
}

### Offer received
AMI$Offer <- sapply(1:N, function(i) which(Z[i, ] == 1))
AMI$Offer <- factor(AMI$Offer)
AMI$Offer <- relevel(AMI$Offer, ref = names(which.max(table(AMI$Offer))))

##### Propensity score #####
### Input variables
Var_XGB <- as.matrix(dummy_cols(AMI[, c('Offer', 'Policy_Type', 'Risk_Group', 'Undershooting_1',
                                        'Undershooting_2', 'Premium_New_Base')], 
                                remove_most_frequent_dummy = FALSE)[, -c(1, 2, 3, 7, 8, 9, 10, 11)])
Var_XGB <- Matrix(Var_XGB, sparse = TRUE)
Labels_XGB <- as.vector(as.numeric(as.character(AMI$Offer))) - 1
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB

### Custom ASAM evaluation function for early stopping
ASAM_feval <- function(preds, dtrain) {
  # Preliminary variables
  N <- nrow(dtrain)
  K <- ncol(dtrain)
  vars <- attr(dtrain, 'data')
  Y_i <- getinfo(dtrain, 'label')
  
  # Determine number of observations in each treatment category
  Z_i <- sparse.model.matrix(~ -1 + factor(Y_i))
  N_t <- apply(Z_i, 2, sum)
  
  # Normalize predictions such that probabilities will sum to one
  preds <- exp(preds) / apply(exp(preds), 1, sum)
  
  # Compute the Average Standardized Absolute Mean (ASAM) difference in the risk factors
  w_mu_t <- t(Z_i / preds) %*% vars / apply(Z_i / preds, 2, sum)
  p_mu <- 1 / N * apply(vars, 2, sum)
  p_sd <- sqrt(diag(t(vars - p_mu) %*% (vars - p_mu)) / (N - 1))
  ASAM <- 1 / T_max * sum(1 / K * apply(t(abs(t(w_mu_t) - p_mu) / p_sd), 1, sum))
  
  # Return the output
  return(list(metric = 'ASAM', value = ASAM))
}

### Settings for Extreme Gradient Boosting (XGBoost)
# General settings
PS_Trees <- 10000
PS_Stop <- 250
seed <- 2019
# Grid search settings
PS_eta <- c(seq(from = 0.01, to = 0.05, by = 0.01), seq(from = 0.1, to = 0.25, by = 0.05), 0.5)
PS_Depth <- c(0, 1, seq(from = 2, to = 10, by = 2), 25, 50)
PS_Min <- c(0, seq(from = 1, to = 5, by = 1), 10, 25, 50) 
PS_row_sample <- c(seq(from = 0.1, to = 1, by = 0.1))
PS_col_sample <- c(seq(from = 0.1, to = 1, by = 0.1))
PS_gamma <- c(0, 0.1, 1, 10, 100)
PS_lambda <- c(0, 0.1, 1, 10, 100)
PS_alpha <- c(0, 0.1, 1, 10, 100)
# Optimal grid
PS_Grid_Opt <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = NA,
  colsample_bytree = NA,
  gamma = NA,
  lambda = NA,
  alpha = NA,
  ASAM = NA,
  Trees = NA
)
# Grid 1 - Learning rate
PS_Grid_1 <- expand.grid(
  eta = PS_eta,
  max_depth = 6,
  min_child_weight = 1,
  subsample = 1,
  colsample_bytree = 1,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  ASAM = NA,
  Trees = NA
)
# Grid 2 - Tree-specific parameters
PS_Grid_2 <- expand.grid(
  eta = NA,
  max_depth = PS_Depth,
  min_child_weight = PS_Min,
  subsample = 1,
  colsample_bytree = 1,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  ASAM = NA,
  Trees = NA
)
# Grid 3 - Stochastic features
PS_Grid_3 <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = PS_row_sample,
  colsample_bytree = PS_col_sample,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  ASAM = NA,
  Trees = NA
)
# Grid 4 - Regularization penalties
PS_Grid_4 <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = NA,
  colsample_bytree = NA,
  gamma = PS_gamma,
  lambda = PS_lambda,
  alpha = PS_alpha,
  ASAM = NA,
  Trees = NA
)

### Grid search 1
for (i in 1:nrow(PS_Grid_1)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(num_class = T_max,
                  objective = 'multi:softprob',
                  booster = 'gbtree', seed = seed,
                  eta = PS_Grid_1$eta[i], 
                  max_depth = PS_Grid_1$max_depth[i], min_child_weight = PS_Grid_1$min_child_weight[i],
                  subsample = PS_Grid_1$subsample[i], colsample_bytree = PS_Grid_1$colsample_bytree[i],
                  gamma = PS_Grid_1$gamma[i], lambda = PS_Grid_1$lambda[i], alpha = PS_Grid_1$alpha[i])
  # Perform cross-validation on the XGBoost model
  PS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                   nfold = 10, nrounds = PS_Trees, early_stopping_rounds = PS_Stop, 
                   verbose = TRUE, maximize = FALSE,
                   feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  PS_Grid_1$ASAM[i] <- min(PS_XGB$evaluation_log$test_ASAM_mean)
  PS_Grid_1$Trees[i] <- which.min(PS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
PS_Grid_Opt$eta <- PS_Grid_1$eta[which.min(PS_Grid_1$ASAM)]
# Adopt these parameters in the next grid searches
PS_Grid_2$eta <- PS_Grid_Opt$eta
PS_Grid_3$eta <- PS_Grid_Opt$eta
PS_Grid_4$eta <- PS_Grid_Opt$eta
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_PS_2 <- which(apply(sapply(1:(ncol(PS_Grid_Opt) - 2), function(j) 
  PS_Grid_2[, j] == PS_Grid_1[which.min(PS_Grid_1$ASAM), j]), 1, sum) == (ncol(PS_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
PS_Grid_2$ASAM[Excl_PS_2] <- PS_Grid_1$ASAM[which.min(PS_Grid_1$ASAM)]
PS_Grid_2$Trees[Excl_PS_2] <- PS_Grid_1$Trees[which.min(PS_Grid_1$ASAM)]

### Grid search 2
for (i in c(1:nrow(PS_Grid_2))[-Excl_PS_2]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(num_class = T_max,
                  objective = 'multi:softprob',
                  booster = 'gbtree', seed = seed,
                  eta = PS_Grid_2$eta[i], 
                  max_depth = PS_Grid_2$max_depth[i], min_child_weight = PS_Grid_2$min_child_weight[i],
                  subsample = PS_Grid_2$subsample[i], colsample_bytree = PS_Grid_2$colsample_bytree[i],
                  gamma = PS_Grid_2$gamma[i], lambda = PS_Grid_2$lambda[i], alpha = PS_Grid_2$alpha[i])
  # Perform cross-validation on the XGBoost model
  PS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                   nfold = 10, nrounds = PS_Trees, early_stopping_rounds = PS_Stop, 
                   verbose = TRUE, maximize = FALSE,
                   feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  PS_Grid_2$ASAM[i] <- min(PS_XGB$evaluation_log$test_ASAM_mean)
  PS_Grid_2$Trees[i] <- which.min(PS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
PS_Grid_Opt$max_depth <- PS_Grid_2$max_depth[which.min(PS_Grid_2$ASAM)]
PS_Grid_Opt$min_child_weight <- PS_Grid_2$min_child_weight[which.min(PS_Grid_2$ASAM)]
# Adopt these parameters in the next grid searches
PS_Grid_3$max_depth <- PS_Grid_Opt$max_depth
PS_Grid_3$min_child_weight <- PS_Grid_Opt$min_child_weight
PS_Grid_4$max_depth <- PS_Grid_Opt$max_depth
PS_Grid_4$min_child_weight <- PS_Grid_Opt$min_child_weight
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_PS_3 <- which(apply(sapply(1:(ncol(PS_Grid_Opt) - 2), function(j) 
  PS_Grid_3[, j] == PS_Grid_2[which.min(PS_Grid_2$ASAM), j]), 1, sum) == (ncol(PS_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
PS_Grid_3$ASAM[Excl_PS_3] <- PS_Grid_2$ASAM[which.min(PS_Grid_2$ASAM)]
PS_Grid_3$Trees[Excl_PS_3] <- PS_Grid_2$Trees[which.min(PS_Grid_2$ASAM)]

### Grid search 3
for (i in c(1:nrow(PS_Grid_3))[-Excl_PS_3]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(num_class = T_max,
                  objective = 'multi:softprob',
                  booster = 'gbtree', seed = seed,
                  eta = PS_Grid_3$eta[i], 
                  max_depth = PS_Grid_3$max_depth[i], min_child_weight = PS_Grid_3$min_child_weight[i],
                  subsample = PS_Grid_3$subsample[i], colsample_bytree = PS_Grid_3$colsample_bytree[i],
                  gamma = PS_Grid_3$gamma[i], lambda = PS_Grid_3$lambda[i], alpha = PS_Grid_3$alpha[i])
  # Perform cross-validation on the XGBoost model
  PS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                   nfold = 10, nrounds = PS_Trees, early_stopping_rounds = PS_Stop, 
                   verbose = TRUE, maximize = FALSE,
                   feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  PS_Grid_3$ASAM[i] <- min(PS_XGB$evaluation_log$test_ASAM_mean)
  PS_Grid_3$Trees[i] <- which.min(PS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
PS_Grid_Opt$subsample <- PS_Grid_3$subsample[which.min(PS_Grid_3$ASAM)]
PS_Grid_Opt$colsample_bytree <- PS_Grid_3$colsample_bytree[which.min(PS_Grid_3$ASAM)]
# Adopt these parameters in the next grid search
PS_Grid_4$subsample <- PS_Grid_Opt$subsample
PS_Grid_4$colsample_bytree <- PS_Grid_Opt$colsample_bytree
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_PS_4 <- which(apply(sapply(1:(ncol(PS_Grid_Opt) - 2), function(j) 
  PS_Grid_4[, j] == PS_Grid_3[which.min(PS_Grid_3$ASAM), j]), 1, sum) == (ncol(PS_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
PS_Grid_4$ASAM[Excl_PS_4] <- PS_Grid_3$ASAM[which.min(PS_Grid_3$ASAM)]
PS_Grid_4$Trees[Excl_PS_4] <- PS_Grid_3$Trees[which.min(PS_Grid_3$ASAM)]

### Grid search 4
for (i in c(1:nrow(PS_Grid_4))[-Excl_PS_4]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(num_class = T_max,
                  objective = 'multi:softprob',
                  booster = 'gbtree', seed = seed,
                  eta = PS_Grid_4$eta[i], 
                  max_depth = PS_Grid_4$max_depth[i], min_child_weight = PS_Grid_4$min_child_weight[i],
                  subsample = PS_Grid_4$subsample[i], colsample_bytree = PS_Grid_4$colsample_bytree[i],
                  gamma = PS_Grid_4$gamma[i], lambda = PS_Grid_4$lambda[i], alpha = PS_Grid_4$alpha[i])
  # Perform cross-validation on the XGBoost model
  PS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                   nfold = 10, nrounds = PS_Trees, early_stopping_rounds = PS_Stop,  
                   verbose = TRUE, maximize = FALSE,
                   feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  PS_Grid_4$ASAM[i] <- min(PS_XGB$evaluation_log$test_ASAM_mean)
  PS_Grid_4$Trees[i] <- which.min(PS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
PS_Grid_Opt$gamma <- PS_Grid_4$gamma[which.min(PS_Grid_4$ASAM)]
PS_Grid_Opt$lambda <- PS_Grid_4$lambda[which.min(PS_Grid_4$ASAM)]
PS_Grid_Opt$alpha <- PS_Grid_4$alpha[which.min(PS_Grid_4$ASAM)]
PS_Grid_Opt$ASAM <- PS_Grid_4$ASAM[which.min(PS_Grid_4$ASAM)]
PS_Grid_Opt$Trees <- PS_Grid_4$Trees[which.min(PS_Grid_4$ASAM)]

### Form the cross-validated XGBoost model
# Optimal number of decision trees
PS_Best <- PS_Grid_Opt$Trees
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Define the parameters of the XGBoost model
par_XGB <- list(num_class = T_max,
                objective = 'multi:softprob', 
                booster = 'gbtree', seed = seed,
                eta = PS_Grid_Opt$eta, 
                max_depth = PS_Grid_Opt$max_depth, min_child_weight = PS_Grid_Opt$min_child_weight,
                subsample = PS_Grid_Opt$subsample, colsample_bytree = PS_Grid_Opt$colsample_bytree,
                gamma = PS_Grid_Opt$gamma, lambda = PS_Grid_Opt$lambda, alpha = PS_Grid_Opt$alpha)
# Estimate the cross-validated XGBoost model
PS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                  nrounds = PS_Best,
                  verbose = TRUE, maximize = FALSE,
                  feval = ASAM_feval)

### Propensity score estimates
# Raw output
Prop_Score <- predict(PS_XGB, newdata = Var_XGB, type = 'response', reshape = TRUE, outputmargin = TRUE)
# Normalize predictions such that probabilities will sum to one
Prop_Score <- exp(Prop_Score) / apply(exp(Prop_Score), 1, sum)

##### Balance #####
### Balance before matching
K <- ncol(Var_XGB)
p_mu <- 1 / N * apply(Var_XGB, 2, sum)
p_sd <- sqrt(diag(t(Var_XGB - p_mu) %*% (Var_XGB - p_mu)) / (N - 1))
mu_t <- as.matrix(t(t(Z) %*% Var_XGB / apply(Z, 2, sum)))
ASAM_t <- abs(mu_t - p_mu) / p_sd
Avg_ASAM_t <- 1 / K * apply(ASAM_t, 2, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score) %*% Var_XGB / apply(Z / Prop_Score, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t <- 1 / K * apply(w_ASAM_t, 2, sum)

### Average propensity scores
# Grid
Step <- 0.40 / 100
T_grid_Adj <- as.numeric(c(-Inf, seq(from = ceiling(Q_bins[1] / Step) * Step,
                                     to = floor(Q_bins[length(Q_bins)]  / Step) * Step,
                                     by = Step), Inf))
T_grid <- as.numeric(c(Q_bins[1],
                       T_grid_Adj[-c(1, length(T_grid_Adj) - 1, length(T_grid_Adj))] + Step / 2, 
                       Q_bins[length(Q_bins)]))
# Aggregation
PS_Dens_Table <- matrix(NA, nrow = length(T_grid_Adj) - 1, ncol = T_max)
for (t in 1:(length(T_grid_Adj) - 1)) {
  # Determine which renewals we should aggregate over
  Table_Indices <- which((AMI$Rate_Change >= T_grid_Adj[t]) & (AMI$Rate_Change < T_grid_Adj[t + 1]))
  # Calculate the average propensity score for these renewals
  if (length(Table_Indices) > 1) {
    PS_Dens_Table[t, ] <- apply(Prop_Score[Table_Indices, ], 2, mean)
  } else {
    PS_Dens_Table[t, ] <- Prop_Score[Table_Indices, ] / length(Table_Indices)
  }
}
# Normalize the average propensity scores such that they sum to one
PS_Dens_Table <- PS_Dens_Table / apply(PS_Dens_Table, 1, sum)

### Store the final results
