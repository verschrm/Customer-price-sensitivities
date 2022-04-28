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

### Generalized propensity score estimates
# Determine number of observations in each discrete treatment category
N_t <- apply(Z, 2, sum)
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_XGB <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_XGB))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_10 <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_10) %*% Var_XGB / apply(Z / GPS_Score_10, 2, sum)))
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

### Generalized propensity score estimates
# Determine number of observations in each discrete treatment category
N_t <- apply(Z, 2, sum)
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_XGB <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_XGB))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_20 <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_20) %*% Var_XGB / apply(Z / GPS_Score_20, 2, sum)))
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

### Generalized propensity score estimates
# Determine number of observations in each discrete treatment category
N_t <- apply(Z, 2, sum)
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_XGB <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_XGB))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_50 <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_50) %*% Var_XGB / apply(Z / GPS_Score_50, 2, sum)))
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

### Generalized propensity score estimates
# Determine number of observations in each discrete treatment category
N_t <- apply(Z, 2, sum)
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_XGB <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_XGB))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_100 <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_100) %*% Var_XGB / apply(Z / GPS_Score_100, 2, sum)))
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

### Generalized propensity score estimates
# Determine number of observations in each discrete treatment category
N_t <- apply(Z, 2, sum)
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_XGB <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_XGB))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_150 <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_150) %*% Var_XGB / apply(Z / GPS_Score_150, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_150 <- 1 / K * apply(w_ASAM_t, 2, sum)

##### 200 categories #####
### Buckets for discrete treatment categories
# 200 buckets
Q_grid <- 1 / 200
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

### Generalized propensity score estimates
# Determine number of observations in each discrete treatment category
N_t <- apply(Z, 2, sum)
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_XGB <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_XGB))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_XGB[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_200 <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_200) %*% Var_XGB / apply(Z / GPS_Score_200, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_200 <- 1 / K * apply(w_ASAM_t, 2, sum)

### Store the final results
