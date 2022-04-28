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
Packages <- c('dplyr', 'fastDummies', 'Matrix', 'lme4', 'gbm', 'xgboost',
              'caret', 'ranger', 'nnet', 'mgcv')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load output from matching procedure
# Load results from Matching_Continuous.R

### Input variables
train_Alt <- AMI[, c('Offer', 'Policy_Type', 'Risk_Group', 'Undershooting_1',
                     'Undershooting_2', 'Premium_New_Base')]
sig2_RF <- pmax(as.vector(1 / (N_t - 1) * t((AMI$Offer - Z_means) ^ 2) %*% Z), 1e-15)
sig2_NN <- pmax(as.vector(1 / (N_t - 1) * t((AMI$Offer - Z_means) ^ 2) %*% Z), 1e-15)

##### Random forest #####
### Custom truncated Normal log-likelihood objective function for random forest
Trunc_Normal_obj_RF <- function(data, lev = NULL, model = NULL) {
  # Preliminary variables
  N <- nrow(data)
  Y_i <- data$obs
  preds <- data$pred
  Z_i <- Z[data$rowIndex, ]
  
  # Determine number of observations in each treatment category
  N_t <- apply(Z_i, 2, sum)
  
  # Adjust heteroskedastic variance estimates
  Var_t <- pmax(as.vector(1 / (N_t - 1) * t((Y_i - preds) ^ 2) %*% Z_i), 1e-15)
  sig2_RF_t <- Var_t
  sig2_t <- as.vector(Z_i %*% as.matrix(Var_t))
  SE <- sqrt(Var_t)
  
  # Pre-calculate some expressions for the objective function
  alpha_i <- sapply(1:T_max, function(t) (Q_bins[t] - preds) / SE[t])
  beta_i <- sapply(1:T_max, function(t) (Q_bins[t + 1] - preds) / SE[t])
  Norm_i <- apply(sapply(1:T_max, function(t) pnorm(beta_i[, t], mean = 0, sd = 1) - 
                           pnorm(alpha_i[, t], mean = 0, sd = 1)), 1, sum)
  
  # Evaluate the objective function
  logloss <- - 1 / N * sum(- log(2 * pi * sqrt(sig2_t)) - ((Y_i - preds) ^ 2) / (2 * sig2_t)  - 
                             log(Norm_i))
  names(logloss) <- 'Trunc.Gaussian'
  
  # Return the output
  return(logloss)
}

### Custom ASAM evaluation function for early stopping for random forest
ASAM_feval_RF <- function(data, lev = NULL, model = NULL) {
  # Preliminary variables
  N <- nrow(data)
  K <- ncol(Var_XGB)
  vars <- Var_XGB[data$rowIndex, ]
  Y_i <- data$obs
  preds <- data$pred
  Z_i <- Z[data$rowIndex, ]
  SE <- sig2_RF
  
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
  names(ASAM) <- 'ASAM'
  
  # Return the output
  return(c(Trunc_Normal_obj_RF(data = data, lev = lev, model = model), ASAM))
}

### Settings for random forest
# General settings
GPS_Trees <- 500
seed <- 2019
# Grid search settings
GPS_max.depth <- c(0, 1, seq(from = 2, to = 10, by = 2), 25, 50)
GPS_min.node.size <- c(0, seq(from = 1, to = 5, by = 1), 10, 25, 50)
GPS_sample.fraction <- c(seq(from = 0.1, to = 1, by = 0.1)) # 1
GPS_mtry <- c(1:(ncol(Var_XGB) - 2)) # 1
GPS_num.tree <- c(1:20, 30, 40, 50, 75, seq(from = 100, to = 5000, by = 50)) # 500
# Optimal grid
GPS_Grid_Opt <- expand.grid(
  max.depth = NA,
  min.node.size = NA,
  sample.fraction = NA,
  mtry = NA,
  num.tree = NA,
  ASAM = NA
)
# Grid 1 - Tree-specific parameters
GPS_Grid_1 <- expand.grid(
  max.depth = GPS_max.depth,
  min.node.size = GPS_min.node.size,
  sample.fraction = 1,
  mtry = 1,
  num.tree = GPS_Trees,
  ASAM = NA
)
# Grid 2 - Stochastic features
GPS_Grid_2 <- expand.grid(
  max.depth = NA,
  min.node.size = NA,
  sample.fraction = GPS_sample.fraction,
  mtry = GPS_mtry,
  num.tree = GPS_Trees,
  ASAM = NA
)
# Grid 3 - Number of trees
GPS_Grid_3 <- expand.grid(
  max.depth = NA,
  min.node.size = NA,
  sample.fraction = NA,
  mtry = NA,
  num.tree = GPS_num.tree,
  ASAM = NA
)

### Grid search 1
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_RF, selectionFunction = 'best')
for (i in 1:nrow(GPS_Grid_1)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  GPS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                    Undershooting_2 + Premium_New_Base, 
                  data = train_Alt, method = 'ranger', metric = 'Trunc.Gaussian',
                  maximize = FALSE, trControl = cv_ctrl,
                  tuneGrid = data.frame(mtry = GPS_Grid_1$mtry[i], splitrule = 'variance',
                                        min.node.size = GPS_Grid_1$min.node.size[i]),
                  max.depth = GPS_Grid_1$max.depth[i],
                  sample.fraction = GPS_Grid_1$sample.fraction[i],
                  num.tree = GPS_Grid_1$num.tree[i], oob.error = FALSE, seed = seed)
  # Store the minimal ASAM obtained
  GPS_Grid_1$ASAM[i] <- GPS_RF$results$ASAM
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$max.depth <- GPS_Grid_1$max.depth[which.min(GPS_Grid_1$ASAM)]
GPS_Grid_Opt$min.node.size <- GPS_Grid_1$min.node.size[which.min(GPS_Grid_1$ASAM)]
# Adopt these parameters in the next grid searches
GPS_Grid_2$max.depth <- GPS_Grid_Opt$max.depth
GPS_Grid_2$min.node.size <- GPS_Grid_Opt$min.node.size
GPS_Grid_3$max.depth <- GPS_Grid_Opt$max.depth
GPS_Grid_3$min.node.size <- GPS_Grid_Opt$min.node.size
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_GPS_2 <- which(apply(sapply(1:(ncol(GPS_Grid_Opt) - 1), function(j) 
  GPS_Grid_2[, j] == GPS_Grid_1[which.min(GPS_Grid_1$ASAM), j]), 1, sum) == (ncol(GPS_Grid_Opt) - 1))
# Store these results for that iteration (as they will return the same outcome)
GPS_Grid_2$ASAM[Excl_GPS_2] <- GPS_Grid_1$ASAM[which.min(GPS_Grid_1$ASAM)]

### Grid search 2
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_RF, selectionFunction = 'best')
for (i in c(1:nrow(GPS_Grid_2))[-Excl_GPS_2]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  GPS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                    Undershooting_2 + Premium_New_Base, 
                  data = train_Alt, method = 'ranger', metric = 'Trunc.Gaussian',
                  maximize = FALSE, trControl = cv_ctrl,
                  tuneGrid = data.frame(mtry = GPS_Grid_2$mtry[i], splitrule = 'variance',
                                        min.node.size = GPS_Grid_2$min.node.size[i]),
                  max.depth = GPS_Grid_2$max.depth[i],
                  sample.fraction = GPS_Grid_2$sample.fraction[i],
                  num.tree = GPS_Grid_2$num.tree[i], oob.error = FALSE, seed = seed)
  # Store the minimal ASAM obtained
  GPS_Grid_2$ASAM[i] <- GPS_RF$results$ASAM
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$sample.fraction <- GPS_Grid_2$sample.fraction[which.min(GPS_Grid_2$ASAM)]
GPS_Grid_Opt$mtry <- GPS_Grid_2$mtry[which.min(GPS_Grid_2$ASAM)]
# Adopt these parameters in the next grid searches
GPS_Grid_3$sample.fraction <- GPS_Grid_Opt$sample.fraction
GPS_Grid_3$mtry <- GPS_Grid_Opt$mtry
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_GPS_3 <- which(apply(sapply(1:(ncol(GPS_Grid_Opt) - 1), function(j) 
  GPS_Grid_3[, j] == GPS_Grid_2[which.min(GPS_Grid_2$ASAM), j]), 1, sum) == (ncol(GPS_Grid_Opt) - 1))
# Store these results for that iteration (as they will return the same outcome)
GPS_Grid_3$ASAM[Excl_GPS_3] <- GPS_Grid_2$ASAM[which.min(GPS_Grid_2$ASAM)]

### Grid search 3
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_RF, selectionFunction = 'best')
for (i in c(1:nrow(GPS_Grid_3))[-Excl_GPS_3]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  GPS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                    Undershooting_2 + Premium_New_Base, 
                  data = train_Alt, method = 'ranger', metric = 'Trunc.Gaussian',
                  maximize = FALSE, trControl = cv_ctrl,
                  tuneGrid = data.frame(mtry = GPS_Grid_3$mtry[i], splitrule = 'variance',
                                        min.node.size = GPS_Grid_3$min.node.size[i]),
                  max.depth = GPS_Grid_3$max.depth[i],
                  sample.fraction = GPS_Grid_3$sample.fraction[i],
                  num.tree = GPS_Grid_3$num.tree[i], oob.error = FALSE, seed = seed)
  # Store the minimal ASAM obtained
  GPS_Grid_3$ASAM[i] <- GPS_RF$results$ASAM
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$num.tree <- GPS_Grid_3$num.tree[which.min(GPS_Grid_3$ASAM)]
GPS_Grid_Opt$ASAM <- GPS_Grid_3$ASAM[which.min(GPS_Grid_3$ASAM)]

### Form the cross-validated random forest model
# Set the cross-validation settings
GPS_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                         summaryFunction = ASAM_feval_RF, selectionFunction = 'best')
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Perform cross-validation on the random forest model
GPS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                  Undershooting_2 + Premium_New_Base, 
                data = train_Alt, method = 'ranger', metric = 'Trunc.Gaussian',
                maximize = FALSE, trControl = GPS_ctrl,
                tuneGrid = data.frame(mtry = GPS_Grid_Opt$mtry, splitrule = 'variance',
                                      min.node.size = GPS_Grid_Opt$min.node.size),
                max.depth = GPS_Grid_Opt$max.depth,
                sample.fraction = GPS_Grid_Opt$sample.fraction,
                num.tree = GPS_Grid_Opt$num.tree, oob.error = FALSE, seed = seed)

### Generalized propensity score estimates
# Retrieve the random forest model predictions
Lin_Pred <- predict(GPS_RF, newdata = train_Alt, type = 'raw')
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_RF <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_RF))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_RF[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_RF[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_RF <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_RF) %*% Var_XGB / apply(Z / GPS_Score_RF, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_RF <- 1 / K * apply(w_ASAM_t, 2, sum)

##### Neural network #####
### Custom truncated Normal log-likelihood objective function for neural network
Trunc_Normal_obj_NN <- function(data, lev = NULL, model = NULL) {
  # Preliminary variables
  N <- nrow(data)
  Y_i <- data$obs
  preds <- data$pred
  Z_i <- Z[data$rowIndex, ]
  
  # Determine number of observations in each treatment category
  N_t <- apply(Z_i, 2, sum)
  
  # Adjust heteroskedastic variance estimates
  Var_t <- pmax(as.vector(1 / (N_t - 1) * t((Y_i - preds) ^ 2) %*% Z_i), 1e-15)
  sig2_NN_t <- Var_t
  sig2_t <- as.vector(Z_i %*% as.matrix(Var_t))
  SE <- sqrt(Var_t)
  
  # Pre-calculate some expressions for the objective function
  alpha_i <- sapply(1:T_max, function(t) (Q_bins[t] - preds) / SE[t])
  beta_i <- sapply(1:T_max, function(t) (Q_bins[t + 1] - preds) / SE[t])
  Norm_i <- apply(sapply(1:T_max, function(t) pnorm(beta_i[, t], mean = 0, sd = 1) - 
                           pnorm(alpha_i[, t], mean = 0, sd = 1)), 1, sum)
  
  # Evaluate the objective function
  logloss <- - 1 / N * sum(- log(2 * pi * sqrt(sig2_t)) - ((Y_i - preds) ^ 2) / (2 * sig2_t)  - 
                             log(Norm_i))
  names(logloss) <- 'Trunc.Gaussian'
  
  # Return the output
  return(logloss)
}

### Custom ASAM evaluation function for early stopping for neural network
ASAM_feval_NN <- function(data, lev = NULL, model = NULL) {
  # Preliminary variables
  N <- nrow(data)
  K <- ncol(Var_XGB)
  vars <- Var_XGB[data$rowIndex, ]
  Y_i <- data$obs
  preds <- data$pred
  Z_i <- Z[data$rowIndex, ]
  SE <- sig2_NN
  
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
  names(ASAM) <- 'ASAM'
  
  # Return the output
  return(c(Trunc_Normal_obj_NN(data = data, lev = lev, model = model), ASAM))
}

### Settings for neural network
# General settings
GPS_Iter <- 100
seed <- 2019
# Grid search settings
GPS_size <- c(1, seq(from = 2, to = 10, by = 2), 25, 50)
GPS_decay <- c(0, 0.1, 1, 10, 100, seq(from = 250, to = 1000, by = 50))
GPS_maxit <- c(1:20, 30, 40, 50, 75, seq(from = 100, to = 5000, by = 50)) # 100
# Optimal grid
GPS_Grid_Opt <- expand.grid(
  size = NA,
  decay = NA,
  maxiter = NA,
  ASAM = NA
)
# Grid 1 - Network-specific parameters
GPS_Grid_1 <- expand.grid(
  size = GPS_size,
  decay = GPS_decay,
  maxiter = GPS_Iter,
  ASAM = NA
)
# Grid 2 - Number of iterations
GPS_Grid_2 <- expand.grid(
  size = NA,
  decay = NA,
  maxiter = GPS_maxit,
  ASAM = NA
)

### Grid search 1
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_NN, selectionFunction = 'best')
for (i in 1:nrow(GPS_Grid_1)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the neural network model
  GPS_NN <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                    Undershooting_2 + Premium_New_Base, 
                  data = train_Alt, method = 'nnet', metric = 'Trunc.Gaussian',
                  preProcess = c('center', 'scale'),
                  maximize = FALSE, trControl = cv_ctrl,
                  tuneGrid = data.frame(size = GPS_Grid_1$size[i],
                                        decay = GPS_Grid_1$decay[i]),
                  trace = FALSE, Hess = FALSE, linout = TRUE, abstol = 0, reltol = 0,
                  maxit = GPS_Grid_1$maxit[i], seed = seed)
  # Store the minimal ASAM obtained
  GPS_Grid_1$ASAM[i] <- GPS_NN$results$ASAM
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$size <- GPS_Grid_1$size[which.min(GPS_Grid_1$ASAM)]
GPS_Grid_Opt$decay <- GPS_Grid_1$decay[which.min(GPS_Grid_1$ASAM)]
# Adopt these parameters in the next grid searches
GPS_Grid_2$size <- GPS_Grid_Opt$size
GPS_Grid_2$decay <- GPS_Grid_Opt$decay
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_GPS_2 <- which(apply(sapply(1:(ncol(GPS_Grid_Opt) - 1), function(j) 
  GPS_Grid_2[, j] == GPS_Grid_1[which.min(GPS_Grid_1$ASAM), j]), 1, sum) == (ncol(GPS_Grid_Opt) - 1))
# Store these results for that iteration (as they will return the same outcome)
GPS_Grid_2$ASAM[Excl_GPS_2] <- GPS_Grid_1$ASAM[which.min(GPS_Grid_1$ASAM)]

### Grid search 2
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_NN, selectionFunction = 'best')
for (i in c(1:nrow(GPS_Grid_2))[-Excl_GPS_2]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the neural network model
  GPS_NN <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                    Undershooting_2 + Premium_New_Base, 
                  data = train_Alt, method = 'nnet', metric = 'Trunc.Gaussian',
                  preProcess = c('center', 'scale'),
                  maximize = FALSE, trControl = cv_ctrl,
                  tuneGrid = data.frame(size = GPS_Grid_2$size[i],
                                        decay = GPS_Grid_2$decay[i]),
                  trace = FALSE, Hess = FALSE, linout = TRUE, abstol = 0, reltol = 0,
                  maxit = GPS_Grid_2$maxit[i], seed = seed)
  # Store the minimal ASAM obtained
  GPS_Grid_2$ASAM[i] <- GPS_RF$results$ASAM
}
# Store the optimal parameters from the grid
GPS_Grid_Opt$maxit <- GPS_Grid_2$maxit[which.min(GPS_Grid_2$ASAM)]
GPS_Grid_Opt$ASAM <- GPS_Grid_2$ASAM[which.min(GPS_Grid_2$ASAM)]

### Form the cross-validated random forest model
# Set the cross-validation settings
GPS_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                         summaryFunction = ASAM_feval_NN, selectionFunction = 'best')
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Perform cross-validation on the random forest model
GPS_NN <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                  Undershooting_2 + Premium_New_Base, 
                data = train_Alt, method = 'nnet', metric = 'Trunc.Gaussian',
                preProcess = c('center', 'scale'),
                maximize = FALSE, trControl = GPS_ctrl,
                tuneGrid = data.frame(size = GPS_Grid_Opt$size,
                                      decay = GPS_Grid_Opt$decay),
                trace = FALSE, Hess = FALSE, linout = TRUE, abstol = 0, reltol = 0,
                maxit = GPS_Grid_Opt$maxit, seed = seed)

### Generalized propensity score estimates
# Retrieve the neural network model predictions
Lin_Pred <- predict(GPS_NN, newdata = train_Alt, type = 'raw')
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_NN <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_NN))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_NN[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_NN[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_NN <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_NN) %*% Var_XGB / apply(Z / GPS_Score_NN, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_NN <- 1 / K * apply(w_ASAM_t, 2, sum)

##### GAM #####
### Form the Generalized Additive Model
GPS_GAM <- gam(formula = Offer ~ s(Undershooting_1, by = Policy_Type, bs = 'cr') + 
                 s(Undershooting_2, by = Policy_Type, bs = 'cr') + 
                 s(Premium_New_Base, by = Policy_Type, bs = 'cr') +
                 s(Undershooting_1, by = Risk_Group, bs = 'cr') + 
                 s(Undershooting_2, by = Risk_Group, bs = 'cr') + 
                 s(Premium_New_Base, by = Risk_Group, bs = 'cr'),
               data = train_Alt, family = gaussian, control = gam.control(maxit = 5000))

### Generalized propensity score estimates
# Retrieve the GAM predictions
Lin_Pred <- predict(GPS_gam_gaus_v8, newdata = train_gaus, type = 'response')
# Determine the corresponding heteroskedastic variance estimates
sd_GPS_GAM <- pmax(as.vector(sqrt(1 / (N_t - 1) * t((AMI$Offer - Lin_Pred) ^ 2) %*% Z)), 1e-15)
sd_GPS_t <- as.vector(Z %*% as.matrix(sd_GPS_GAM))
# Aggregate and normalize the predictions over the treatment intervals
GPS_Score <- sapply(1:T_max, function(t) 
  pnorm((Q_bins[t + 1] - Lin_Pred) / sd_GPS_GAM[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Lin_Pred) / sd_GPS_GAM[t], mean = 0, sd = 1))
GPS_Norm <- apply(GPS_Score, 1, sum)
GPS_Score_GAM <- GPS_Score / GPS_Norm

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_GAM) %*% Var_XGB / apply(Z / GPS_Score_GAM, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_GAM <- 1 / K * apply(w_ASAM_t, 2, sum)

### Store the final results
