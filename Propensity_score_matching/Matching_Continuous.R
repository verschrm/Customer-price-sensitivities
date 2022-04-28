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

##### Generalized propensity score #####
### Input variables
Var_XGB <- as.matrix(dummy_cols(AMI[, c('Offer', 'Policy_Type', 'Risk_Group', 'Undershooting_1',
                                        'Undershooting_2', 'Premium_New_Base')], 
                                remove_most_frequent_dummy = FALSE)[, -c(1, 2, 3)])
Var_XGB <- Matrix(Var_XGB, sparse = TRUE)
K <- ncol(Var_XGB)
Labels_XGB <- as.vector(AMI$Offer)
# Determine number of observations in each discrete treatment category
N_t <- apply(Z, 2, sum)
# Determine the initial points for heteroskedastic variance estimates
Z_means <- as.vector(sapply(1:T_max, function(t) (Q_bins[t] + Q_bins[t + 1]) / 2) * Z)
Z_means <- Z_means[which(Z_means != 0)]
sig2_XGB <- pmax(as.vector(1 / (N_t - 1) * t((AMI$Offer - Z_means) ^ 2) %*% Z), 1e-15)
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB
attr(dtrain_XGB, 'Z_i') <- Z
attr(dtrain_XGB, 'Var_t') <- sig2_XGB
attr(dtrain_XGB, 'Var_init') <- sig2_XGB

### Adjusted buckets
Q_bins_Adj <- Q_bins
Q_bins_Adj[1] <- -Inf
Q_bins_Adj[length(Q_bins_Adj)] <- Inf

### Custom truncated Normal log-likelihood objective function
Trunc_Normal_obj <- function(preds, dtrain) {
  # Preliminary variables
  Y_i <- getinfo(dtrain, 'label')
  N <- length(Y_i)
  Z_i <- attr(dtrain, 'Z_i')
  
  # Determine number of observations in each treatment category
  N_t <- apply(Z_i, 2, sum)
  
  # Adjust heteroskedastic variance estimates
  Var_t <- as.vector(1 / (N_t - 1) * t((Y_i - preds) ^ 2) %*% Z_i)
  attr(dtrain_XGB, 'Var_t') <- Var_t
  sig2_t <- as.vector(Z_i %*% as.matrix(Var_t))
  SE <- sqrt(Var_t)
  
  # Pre-calculate some expressions for the Gradient and Hessian of the objective function
  alpha_i <- sapply(1:T_max, function(t) (Q_bins[t] - preds) / SE[t])
  beta_i <- sapply(1:T_max, function(t) (Q_bins[t + 1] - preds) / SE[t])
  Norm_i <- apply(sapply(1:T_max, function(t) pnorm(beta_i[, t], mean = 0, sd = 1) - 
                           pnorm(alpha_i[, t], mean = 0, sd = 1)), 1, sum)
  phi_alpha_i <- sapply(1:T_max, function(t) dnorm(alpha_i[, t], mean = 0, sd = 1))
  phi_beta_i <- sapply(1:T_max, function(t) dnorm(beta_i[, t], mean = 0, sd = 1))
  
  # Combine previous calculations into the Gradient and Hessian of the objective function
  Gradient <- - 1 / N * ((Y_i - preds) / sig2_t - apply(phi_alpha_i - phi_beta_i, 1, sum) / Norm_i)
  Hessian <- - 1 / N *(-1 / sig2_t - apply(alpha_i * phi_alpha_i - beta_i * phi_beta_i, 1, function(i) 
    sum(i / SE)) / Norm_i + apply(phi_alpha_i - phi_beta_i, 1, sum))
  
  # Return the output
  return(list(grad = as.vector(Gradient), hess = as.vector(Hessian)))
}

### Custom ASAM evaluation function for early stopping
ASAM_feval <- function(preds, dtrain) {
  # Preliminary variables
  N <- nrow(dtrain)
  K <- ncol(dtrain)
  vars <- attr(dtrain, 'data')
  Y_i <- getinfo(dtrain, 'label')
  Z_i <- attr(dtrain, 'Z_i')
  SE <- sqrt(attr(dtrain_XGB, 'Var_t'))
  
  # Determine number of observations in each treatment category
  N_t <- apply(Z_i, 2, sum)
  
  # Aggregate and normalize the predictions over the treatment intervals
  preds <- sapply(1:T_max, function(t) pnorm((Q_bins[t + 1] - preds) / SE[t], mean = 0, sd = 1) - 
                    pnorm((Q_bins[t] - preds) / SE[t], mean = 0, sd = 1))
  preds <- preds / apply(preds, 1, sum)
  
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
GPS_Trees <- 10000
GPS_Stop <- 250
seed <- 2019
# Grid search settings
GPS_eta <- c(seq(from = 0.01, to = 0.05, by = 0.01), seq(from = 0.1, to = 0.25, by = 0.05), 0.5)
GPS_Depth <- c(0, 1, seq(from = 2, to = 10, by = 2), 25, 50) # 6
GPS_Min <- c(0, seq(from = 1, to = 5, by = 1), 10, 25, 50) # 1
GPS_row_sample <- c(seq(from = 0.1, to = 1, by = 0.1)) # 1
GPS_col_sample <- c(seq(from = 0.1, to = 1, by = 0.1)) # 1
GPS_gamma <- c(0, 0.1, 1, 10, 100) # 0
GPS_lambda <- c(0, 0.1, 1, 10, 100) # 1
GPS_alpha <- c(0, 0.1, 1, 10, 100) # 0
# Optimal grid
GPS_Grid_Opt <- expand.grid(
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
GPS_Grid_1 <- expand.grid(
  eta = GPS_eta,
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
GPS_Grid_2 <- expand.grid(
  eta = NA,
  max_depth = GPS_Depth,
  min_child_weight = GPS_Min,
  subsample = 1,
  colsample_bytree = 1,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  ASAM = NA,
  Trees = NA
)
# Grid 3 - Stochastic features
GPS_Grid_3 <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = GPS_row_sample,
  colsample_bytree = GPS_col_sample,
  gamma = 0,
  lambda = 1,
  alpha = 0,
  ASAM = NA,
  Trees = NA
)
# Grid 4 - Regularization penalties
GPS_Grid_4 <- expand.grid(
  eta = NA,
  max_depth = NA,
  min_child_weight = NA,
  subsample = NA,
  colsample_bytree = NA,
  gamma = GPS_gamma,
  lambda = GPS_lambda,
  alpha = GPS_alpha,
  ASAM = NA,
  Trees = NA
)

### Grid search 1
for (i in 1:nrow(GPS_Grid_1)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(booster = 'gbtree', seed = seed, eta = GPS_Grid_1$eta[i], 
                  max_depth = GPS_Grid_1$max_depth[i], 
                  min_child_weight = GPS_Grid_1$min_child_weight[i],
                  subsample = GPS_Grid_1$subsample[i], 
                  colsample_bytree = GPS_Grid_1$colsample_bytree[i],
                  gamma = GPS_Grid_1$gamma[i], lambda = GPS_Grid_1$lambda[i], 
                  alpha = GPS_Grid_1$alpha[i])
  # Perform cross-validation on the XGBoost model
  GPS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                    nfold = 10, nrounds = GPS_Trees, early_stopping_rounds = GPS_Stop, 
                    verbose = TRUE, maximize = FALSE,
                    obj = Trunc_Normal_obj, base_score = (Q_bins[1] + Q_bins[T_max + 1]) / 2,
                    feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  GPS_Grid_1$ASAM[i] <- min(GPS_XGB$evaluation_log$test_ASAM_mean)
  GPS_Grid_1$Trees[i] <- which.min(GPS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$eta <- GPS_Grid_1$eta[which.min(GPS_Grid_1$ASAM)]
# Adopt these parameters in the next grid searches
GPS_Grid_2$eta <- GPS_Grid_Opt$eta
GPS_Grid_3$eta <- GPS_Grid_Opt$eta
GPS_Grid_4$eta <- GPS_Grid_Opt$eta
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_GPS_2 <- which(apply(sapply(1:(ncol(GPS_Grid_Opt) - 2), function(j) 
  GPS_Grid_2[, j] == GPS_Grid_1[which.min(GPS_Grid_1$ASAM), j]), 1, sum) == (ncol(GPS_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
GPS_Grid_2$ASAM[Excl_GPS_2] <- GPS_Grid_1$ASAM[which.min(GPS_Grid_1$ASAM)]
GPS_Grid_2$Trees[Excl_GPS_2] <- GPS_Grid_1$Trees[which.min(GPS_Grid_1$ASAM)]

### Grid search 2
for (i in c(1:nrow(GPS_Grid_2))[-Excl_GPS_2]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(booster = 'gbtree', seed = seed, eta = GPS_Grid_2$eta[i], 
                  max_depth = GPS_Grid_2$max_depth[i], 
                  min_child_weight = GPS_Grid_2$min_child_weight[i],
                  subsample = GPS_Grid_2$subsample[i], 
                  colsample_bytree = GPS_Grid_2$colsample_bytree[i],
                  gamma = GPS_Grid_2$gamma[i], lambda = GPS_Grid_2$lambda[i], 
                  alpha = GPS_Grid_2$alpha[i])
  # Perform cross-validation on the XGBoost model
  GPS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                    nfold = 10, nrounds = GPS_Trees, early_stopping_rounds = GPS_Stop, 
                    verbose = TRUE, maximize = FALSE,
                    obj = Trunc_Normal_obj, base_score = (Q_bins[1] + Q_bins[T_max + 1]) / 2,
                    feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  GPS_Grid_2$ASAM[i] <- min(GPS_XGB$evaluation_log$test_ASAM_mean)
  GPS_Grid_2$Trees[i] <- which.min(GPS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$max_depth <- GPS_Grid_2$max_depth[which.min(GPS_Grid_2$ASAM)]
GPS_Grid_Opt$min_child_weight <- GPS_Grid_2$min_child_weight[which.min(GPS_Grid_2$ASAM)]
# Adopt these parameters in the next grid searches
GPS_Grid_3$max_depth <- GPS_Grid_Opt$max_depth
GPS_Grid_3$min_child_weight <- GPS_Grid_Opt$min_child_weight
GPS_Grid_4$max_depth <- GPS_Grid_Opt$max_depth
GPS_Grid_4$min_child_weight <- GPS_Grid_Opt$min_child_weight
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_GPS_3 <- which(apply(sapply(1:(ncol(GPS_Grid_Opt) - 2), function(j) 
  GPS_Grid_3[, j] == GPS_Grid_2[which.min(GPS_Grid_2$ASAM), j]), 1, sum) == (ncol(GPS_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
GPS_Grid_3$ASAM[Excl_GPS_3] <- GPS_Grid_2$ASAM[which.min(GPS_Grid_2$ASAM)]
GPS_Grid_3$Trees[Excl_GPS_3] <- GPS_Grid_2$Trees[which.min(GPS_Grid_2$ASAM)]

### Grid search 3
for (i in c(1:nrow(GPS_Grid_3))[-Excl_GPS_3]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(booster = 'gbtree', seed = seed, eta = GPS_Grid_3$eta[i], 
                  max_depth = GPS_Grid_3$max_depth[i], 
                  min_child_weight = GPS_Grid_3$min_child_weight[i],
                  subsample = GPS_Grid_3$subsample[i], 
                  colsample_bytree = GPS_Grid_3$colsample_bytree[i],
                  gamma = GPS_Grid_3$gamma[i], lambda = GPS_Grid_3$lambda[i], 
                  alpha = GPS_Grid_3$alpha[i])
  # Perform cross-validation on the XGBoost model
  GPS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                    nfold = 10, nrounds = GPS_Trees, early_stopping_rounds = GPS_Stop, 
                    verbose = TRUE, maximize = FALSE,
                    obj = Trunc_Normal_obj, base_score = (Q_bins[1] + Q_bins[T_max + 1]) / 2,
                    feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  GPS_Grid_3$ASAM[i] <- min(GPS_XGB$evaluation_log$test_ASAM_mean)
  GPS_Grid_3$Trees[i] <- which.min(GPS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$subsample <- GPS_Grid_3$subsample[which.min(GPS_Grid_3$ASAM)]
GPS_Grid_Opt$colsample_bytree <- GPS_Grid_3$colsample_bytree[which.min(GPS_Grid_3$ASAM)]
# Adopt these parameters in the next grid search
GPS_Grid_4$subsample <- GPS_Grid_Opt$subsample
GPS_Grid_4$colsample_bytree <- GPS_Grid_Opt$colsample_bytree
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_GPS_4 <- which(apply(sapply(1:(ncol(GPS_Grid_Opt) - 2), function(j) 
  GPS_Grid_4[, j] == GPS_Grid_3[which.min(GPS_Grid_3$ASAM), j]), 1, sum) == (ncol(GPS_Grid_Opt) - 2))
# Store these results for that iteration (as they will return the same outcome)
GPS_Grid_4$ASAM[Excl_GPS_4] <- GPS_Grid_3$ASAM[which.min(GPS_Grid_3$ASAM)]
GPS_Grid_4$Trees[Excl_GPS_4] <- GPS_Grid_3$Trees[which.min(GPS_Grid_3$ASAM)]

### Grid search 4
for (i in c(1:nrow(GPS_Grid_4))[-Excl_GPS_4]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Define the parameters of the XGBoost model
  par_XGB <- list(booster = 'gbtree', seed = seed, eta = GPS_Grid_4$eta[i], 
                  max_depth = GPS_Grid_4$max_depth[i], 
                  min_child_weight = GPS_Grid_4$min_child_weight[i],
                  subsample = GPS_Grid_4$subsample[i], 
                  colsample_bytree = GPS_Grid_4$colsample_bytree[i],
                  gamma = GPS_Grid_4$gamma[i], lambda = GPS_Grid_4$lambda[i], 
                  alpha = GPS_Grid_4$alpha[i])
  # Perform cross-validation on the XGBoost model
  GPS_XGB <- xgb.cv(data = dtrain_XGB, params = par_XGB,
                    nfold = 10, nrounds = GPS_Trees, early_stopping_rounds = GPS_Stop,  
                    verbose = TRUE, maximize = FALSE,
                    obj = Trunc_Normal_obj, base_score = (Q_bins[1] + Q_bins[T_max + 1]) / 2,
                    feval = ASAM_feval)
  # Store the minimal ASAM obtained and the corresponding number of decision trees used
  GPS_Grid_4$ASAM[i] <- min(GPS_XGB$evaluation_log$test_ASAM_mean)
  GPS_Grid_4$Trees[i] <- which.min(GPS_XGB$evaluation_log$test_ASAM_mean)
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$gamma <- GPS_Grid_4$gamma[which.min(GPS_Grid_4$ASAM)]
GPS_Grid_Opt$lambda <- GPS_Grid_4$lambda[which.min(GPS_Grid_4$ASAM)]
GPS_Grid_Opt$alpha <- GPS_Grid_4$alpha[which.min(GPS_Grid_4$ASAM)]
GPS_Grid_Opt$ASAM <- GPS_Grid_4$ASAM[which.min(GPS_Grid_4$ASAM)]
GPS_Grid_Opt$Trees <- GPS_Grid_4$Trees[which.min(GPS_Grid_4$ASAM)]

### Form the cross-validated XGBoost model
# Optimal number of decision trees
GPS_Best <- GPS_Grid_Opt$Trees
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Define the parameters of the XGBoost model
par_XGB <- list(booster = 'gbtree', seed = seed, eta = GPS_Grid_Opt$eta, 
                max_depth = GPS_Grid_Opt$max_depth, 
                min_child_weight = GPS_Grid_Opt$min_child_weight,
                subsample = GPS_Grid_Opt$subsample, 
                colsample_bytree = GPS_Grid_Opt$colsample_bytree,
                gamma = GPS_Grid_Opt$gamma, lambda = GPS_Grid_Opt$lambda, 
                alpha = GPS_Grid_Opt$alpha)
# Estimate the cross-validated XGBoost model
GPS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                   nrounds = GPS_Best,
                   verbose = TRUE, maximize = FALSE,
                   obj = Trunc_Normal_obj, base_score = (Q_bins[1] + Q_bins[T_max + 1]) / 2,
                   feval = ASAM_feval)

### Generalized propensity score estimates
# Retrieve the XGBoost model predictions
Lin_Pred <- predict(GPS_XGB, newdata = as.matrix(Var_XGB), type = 'response')
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_XGB <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_XGB))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score <- GPS_Score / GPS_Norm
# Determine the observed generalized propensity scores
AMI$GPS <- dnorm((AMI$Offer - Lin_Pred) / sd_GPS_t, mean = 0, sd = 1) / 
  GPS_Norm / sd_GPS_t

##### Balance #####
### Balance before matching
K <- ncol(Var_XGB)
p_mu <- 1 / N * apply(Var_XGB, 2, sum)
p_sd <- sqrt(diag(t(Var_XGB - p_mu) %*% (Var_XGB - p_mu)) / (N - 1))
mu_t <- as.matrix(t(t(Z) %*% Var_XGB / apply(Z, 2, sum)))
ASAM_t <- abs(mu_t - p_mu) / p_sd
Avg_ASAM_t <- 1 / K * apply(ASAM_t, 2, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score) %*% Var_XGB / apply(Z / GPS_Score, 2, sum)))
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
GPS_Dens_Table <- matrix(NA, nrow = length(T_grid_Adj) - 1, ncol = T_max)
for (t in 1:(length(T_grid_Adj) - 1)) {
  # Determine which renewals we should aggregate over
  Table_Indices <- which((AMI$Rate_Change >= T_grid_Adj[t]) & (AMI$Rate_Change < T_grid_Adj[t + 1]))
  # Calculate the average propensity score for these renewals
  if (length(Table_Indices) > 1) {
    GPS_Dens_Table[t, ] <- apply(GPS_Score[Table_Indices, ], 2, mean)
  } else {
    GPS_Dens_Table[t, ] <- GPS_Score[Table_Indices, ] / length(Table_Indices)
  }
}
# Normalize the average propensity scores such that they sum to one
GPS_Dens_Table <- GPS_Dens_Table / apply(GPS_Dens_Table, 1, sum)

### Store the final results
