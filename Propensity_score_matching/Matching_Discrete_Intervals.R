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
# Load results from Matching_Discrete.R

##### 10 categories #####
### Buckets for discrete treatment categories
# 10 buckets
Q_grid <- 1 / 10
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

### Input variables
Labels_XGB <- as.vector(as.numeric(as.character(AMI$Offer))) - 1
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB

### Optimized XGBoost model
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
# Estimate the optimized XGBoost model
PS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                  nrounds = PS_Best,
                  verbose = TRUE, maximize = FALSE,
                  feval = ASAM_feval)

### Propensity score estimates
# Raw output
Prop_Score <- predict(PS_XGB, newdata = Var_XGB, type = 'response', reshape = TRUE, outputmargin = TRUE)
# Normalize predictions such that probabilities will sum to one
Prop_Score_10 <- exp(Prop_Score) / apply(exp(Prop_Score), 1, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_10) %*% Var_XGB / apply(Z / Prop_Score_10, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_10 <- 1 / K * apply(w_ASAM_t, 2, sum)

##### 20 categories #####
### Buckets for discrete treatment categories
# 20 buckets
Q_grid <- 1 / 20
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

### Input variables
Labels_XGB <- as.vector(as.numeric(as.character(AMI$Offer))) - 1
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB

### Optimized XGBoost model
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
# Estimate the optimized XGBoost model
PS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                  nrounds = PS_Best,
                  verbose = TRUE, maximize = FALSE,
                  feval = ASAM_feval)

### Propensity score estimates
# Raw output
Prop_Score <- predict(PS_XGB, newdata = Var_XGB, type = 'response', reshape = TRUE, outputmargin = TRUE)
# Normalize predictions such that probabilities will sum to one
Prop_Score_20 <- exp(Prop_Score) / apply(exp(Prop_Score), 1, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_20) %*% Var_XGB / apply(Z / Prop_Score_20, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_20 <- 1 / K * apply(w_ASAM_t, 2, sum)

##### 50 categories #####
### Buckets for discrete treatment categories
# 50 buckets
Q_grid <- 1 / 50
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

### Input variables
Labels_XGB <- as.vector(as.numeric(as.character(AMI$Offer))) - 1
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB

### Optimized XGBoost model
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
# Estimate the optimized XGBoost model
PS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                  nrounds = PS_Best,
                  verbose = TRUE, maximize = FALSE,
                  feval = ASAM_feval)

### Propensity score estimates
# Raw output
Prop_Score <- predict(PS_XGB, newdata = Var_XGB, type = 'response', reshape = TRUE, outputmargin = TRUE)
# Normalize predictions such that probabilities will sum to one
Prop_Score_50 <- exp(Prop_Score) / apply(exp(Prop_Score), 1, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_50) %*% Var_XGB / apply(Z / Prop_Score_50, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_50 <- 1 / K * apply(w_ASAM_t, 2, sum)

##### 100 categories #####
### Buckets for discrete treatment categories
# 100 buckets
Q_grid <- 1 / 100
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

### Input variables
Labels_XGB <- as.vector(as.numeric(as.character(AMI$Offer))) - 1
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB

### Optimized XGBoost model
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
# Estimate the optimized XGBoost model
PS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                  nrounds = PS_Best,
                  verbose = TRUE, maximize = FALSE,
                  feval = ASAM_feval)

### Propensity score estimates
# Raw output
Prop_Score <- predict(PS_XGB, newdata = Var_XGB, type = 'response', reshape = TRUE, outputmargin = TRUE)
# Normalize predictions such that probabilities will sum to one
Prop_Score_100 <- exp(Prop_Score) / apply(exp(Prop_Score), 1, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_100) %*% Var_XGB / apply(Z / Prop_Score_100, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_100 <- 1 / K * apply(w_ASAM_t, 2, sum)

##### 150 categories #####
### Buckets for discrete treatment categories
# 150 buckets
Q_grid <- 1 / 150
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

### Input variables
Labels_XGB <- as.vector(as.numeric(as.character(AMI$Offer))) - 1
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB

### Optimized XGBoost model
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
# Estimate the optimized XGBoost model
PS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                  nrounds = PS_Best,
                  verbose = TRUE, maximize = FALSE,
                  feval = ASAM_feval)

### Propensity score estimates
# Raw output
Prop_Score <- predict(PS_XGB, newdata = Var_XGB, type = 'response', reshape = TRUE, outputmargin = TRUE)
# Normalize predictions such that probabilities will sum to one
Prop_Score_150 <- exp(Prop_Score) / apply(exp(Prop_Score), 1, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_150) %*% Var_XGB / apply(Z / Prop_Score_150, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_150 <- 1 / K * apply(w_ASAM_t, 2, sum)

##### 200 categories #####
### Buckets for discrete treatment categories
# 200 buckets
Q_grid <- 1 / 10
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

### Input variables
Labels_XGB <- as.vector(as.numeric(as.character(AMI$Offer))) - 1
# Assign to matrix object for XGBoost
dtrain_XGB <- xgb.DMatrix(data = as.matrix(Var_XGB), label = Labels_XGB)
attr(dtrain_XGB, 'data') <- Var_XGB

### Optimized XGBoost model
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
# Estimate the optimized XGBoost model
PS_XGB <- xgboost(data = dtrain_XGB, params = par_XGB, 
                  nrounds = PS_Best,
                  verbose = TRUE, maximize = FALSE,
                  feval = ASAM_feval)

### Propensity score estimates
# Raw output
Prop_Score <- predict(PS_XGB, newdata = Var_XGB, type = 'response', reshape = TRUE, outputmargin = TRUE)
# Normalize predictions such that probabilities will sum to one
Prop_Score_200 <- exp(Prop_Score) / apply(exp(Prop_Score), 1, sum)

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_200) %*% Var_XGB / apply(Z / Prop_Score_200, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_200 <- 1 / K * apply(w_ASAM_t, 2, sum)

### Store the final results
