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
# Load results from Matching_Discrete.R

### Input variables
train_Alt <- AMI[, c('Offer', 'Policy_Type', 'Risk_Group', 'Undershooting_1',
                     'Undershooting_2', 'Premium_New_Base')]

##### Random forest #####
### Custom ASAM evaluation function for early stopping for random forest
ASAM_feval_RF <- function(data, lev = NULL, model = NULL) {
  # Preliminary variables
  N <- nrow(data)
  K <- ncol(Var_XGB)
  vars <- Var_XGB[data$rowIndex, ]
  Y_i <- data$obs
  preds <- as.matrix(data[, c(lev)])
  
  # Determine number of observations in each treatment category
  Z_i <- sparse.model.matrix(~ -1 + factor(Y_i))
  N_t <- apply(Z_i, 2, sum)
  
  # Compute the Average Standardized Absolute Mean (ASAM) difference in the risk factors
  w_mu_t <- t(Z_i / preds) %*% vars / apply(Z_i / preds, 2, sum)
  p_mu <- 1 / N * apply(vars, 2, sum)
  p_sd <- sqrt(diag(t(vars - p_mu) %*% (vars - p_mu)) / (N - 1))
  ASAM <- 1 / T_max * sum(1 / K * apply(t(abs(t(w_mu_t) - p_mu) / p_sd), 1, sum))
  names(ASAM) <- 'ASAM'
  
  # Return the output
  return(c(mnLogLoss(data = data, lev = lev, model = model), ASAM))
}

### Settings for random forest
# General settings
PS_Trees <- 500
seed <- 2019
# Grid search settings
PS_max.depth <- c(0, 1, seq(from = 2, to = 10, by = 2), 25, 50)
PS_min.node.size <- c(0, seq(from = 1, to = 5, by = 1), 10, 25, 50)
PS_sample.fraction <- c(seq(from = 0.1, to = 1, by = 0.1)) # 1
PS_mtry <- c(1:(ncol(Var_XGB) - 2)) # 1
PS_num.tree <- c(1:20, 30, 40, 50, 75, seq(from = 100, to = 5000, by = 50)) # 500
# Optimal grid
PS_Grid_Opt <- expand.grid(
  max.depth = NA,
  min.node.size = NA,
  sample.fraction = NA,
  mtry = NA,
  num.tree = NA,
  ASAM = NA
)
# Grid 1 - Tree-specific parameters
PS_Grid_1 <- expand.grid(
  max.depth = PS_max.depth,
  min.node.size = PS_min.node.size,
  sample.fraction = 1,
  mtry = 1,
  num.tree = PS_Trees,
  ASAM = NA
)
# Grid 2 - Stochastic features
PS_Grid_2 <- expand.grid(
  max.depth = NA,
  min.node.size = NA,
  sample.fraction = PS_sample.fraction,
  mtry = PS_mtry,
  num.tree = PS_Trees,
  ASAM = NA
)
# Grid 3 - Number of trees
PS_Grid_3 <- expand.grid(
  max.depth = NA,
  min.node.size = NA,
  sample.fraction = NA,
  mtry = NA,
  num.tree = PS_num.tree,
  ASAM = NA
)

### Grid search 1
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_RF, selectionFunction = 'best',
                        classProbs = TRUE)
for (i in 1:nrow(PS_Grid_1)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  PS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                   Undershooting_2 + Premium_New_Base, 
                 data = train_Alt, method = 'ranger', metric = 'logLoss',
                 maximize = FALSE, trControl = cv_ctrl,
                 tuneGrid = data.frame(mtry = PS_Grid_1$mtry[i], splitrule = 'gini',
                                       min.node.size = PS_Grid_1$min.node.size[i]),
                 max.depth = PS_Grid_1$max.depth[i],
                 sample.fraction = PS_Grid_1$sample.fraction[i],
                 num.tree = PS_Grid_1$num.tree[i], oob.error = FALSE, seed = seed)
  # Store the minimal ASAM obtained
  PS_Grid_1$ASAM[i] <- PS_RF$results$ASAM
}
# Store the optimal parameters from the grid
PS_Grid_Opt$max.depth <- PS_Grid_1$max.depth[which.min(PS_Grid_1$ASAM)]
PS_Grid_Opt$min.node.size <- PS_Grid_1$min.node.size[which.min(PS_Grid_1$ASAM)]
# Adopt these parameters in the next grid searches
PS_Grid_2$max.depth <- PS_Grid_Opt$max.depth
PS_Grid_2$min.node.size <- PS_Grid_Opt$min.node.size
PS_Grid_3$max.depth <- PS_Grid_Opt$max.depth
PS_Grid_3$min.node.size <- PS_Grid_Opt$min.node.size
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_PS_2 <- which(apply(sapply(1:(ncol(PS_Grid_Opt) - 1), function(j) 
  PS_Grid_2[, j] == PS_Grid_1[which.min(PS_Grid_1$ASAM), j]), 1, sum) == (ncol(PS_Grid_Opt) - 1))
# Store these results for that iteration (as they will return the same outcome)
PS_Grid_2$ASAM[Excl_PS_2] <- PS_Grid_1$ASAM[which.min(PS_Grid_1$ASAM)]

### Grid search 2
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_RF, selectionFunction = 'best',
                        classProbs = TRUE)
for (i in c(1:nrow(PS_Grid_2))[-Excl_PS_2]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  PS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                   Undershooting_2 + Premium_New_Base, 
                 data = train_Alt, method = 'ranger', metric = 'logLoss',
                 maximize = FALSE, trControl = cv_ctrl,
                 tuneGrid = data.frame(mtry = PS_Grid_2$mtry[i], splitrule = 'gini',
                                       min.node.size = PS_Grid_2$min.node.size[i]),
                 max.depth = PS_Grid_2$max.depth[i],
                 sample.fraction = PS_Grid_2$sample.fraction[i],
                 num.tree = PS_Grid_2$num.tree[i], oob.error = FALSE, seed = seed)
  # Store the minimal ASAM obtained
  PS_Grid_2$ASAM[i] <- PS_RF$results$ASAM
}
# Store the optimal parameters from the grid
PS_Grid_Opt$sample.fraction <- PS_Grid_2$sample.fraction[which.min(PS_Grid_2$ASAM)]
PS_Grid_Opt$mtry <- PS_Grid_2$mtry[which.min(PS_Grid_2$ASAM)]
# Adopt these parameters in the next grid searches
PS_Grid_3$sample.fraction <- PS_Grid_Opt$sample.fraction
PS_Grid_3$mtry <- PS_Grid_Opt$mtry
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_PS_3 <- which(apply(sapply(1:(ncol(PS_Grid_Opt) - 1), function(j) 
  PS_Grid_3[, j] == PS_Grid_2[which.min(PS_Grid_2$ASAM), j]), 1, sum) == (ncol(PS_Grid_Opt) - 1))
# Store these results for that iteration (as they will return the same outcome)
PS_Grid_3$ASAM[Excl_PS_3] <- PS_Grid_2$ASAM[which.min(PS_Grid_2$ASAM)]

### Grid search 3
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_RF, selectionFunction = 'best',
                        classProbs = TRUE)
for (i in c(1:nrow(PS_Grid_3))[-Excl_PS_3]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  PS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                   Undershooting_2 + Premium_New_Base, 
                 data = train_Alt, method = 'ranger', metric = 'logLoss',
                 maximize = FALSE, trControl = cv_ctrl,
                 tuneGrid = data.frame(mtry = PS_Grid_3$mtry[i], splitrule = 'gini',
                                       min.node.size = PS_Grid_3$min.node.size[i]),
                 max.depth = PS_Grid_3$max.depth[i],
                 sample.fraction = PS_Grid_3$sample.fraction[i],
                 num.tree = PS_Grid_3$num.tree[i], oob.error = FALSE, seed = seed)
  # Store the minimal ASAM obtained
  PS_Grid_3$ASAM[i] <- PS_RF$results$ASAM
}
# Store the optimal parameters from the grid
PS_Grid_Opt$num.tree <- PS_Grid_3$num.tree[which.min(PS_Grid_3$ASAM)]
PS_Grid_Opt$ASAM <- PS_Grid_3$ASAM[which.min(PS_Grid_3$ASAM)]

### Form the cross-validated random forest model
# Set the cross-validation settings
PS_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_RF, selectionFunction = 'best',
                        classProbs = TRUE)
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Perform cross-validation on the random forest model
PS_RF <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                 Undershooting_2 + Premium_New_Base, 
               data = train_Alt, method = 'ranger', metric = 'logLoss',
               maximize = FALSE, trControl = PS_ctrl,
               tuneGrid = data.frame(mtry = PS_Grid_Opt$mtry, splitrule = 'gini',
                                     min.node.size = PS_Grid_Opt$min.node.size),
               max.depth = PS_Grid_Opt$max.depth,
               sample.fraction = PS_Grid_Opt$sample.fraction,
               num.tree = PS_Grid_Opt$num.tree, oob.error = FALSE, seed = seed)

### Propensity score estimates
Prop_Score_RF <- predict(PS_RF, newdata = train_Alt, type = 'prob')

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_RF) %*% Var_XGB / apply(Z / Prop_Score_RF, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_RF <- 1 / K * apply(w_ASAM_t, 2, sum)

##### Neural network #####
### Custom ASAM evaluation function for early stopping for neural network
ASAM_feval_NN <- function(data, lev = NULL, model = NULL) {
  # Preliminary variables
  N <- nrow(data)
  K <- ncol(Var_XGB)
  vars <- Var_XGB[data$rowIndex, ]
  Y_i <- data$obs
  preds <- as.matrix(data[, c(lev)])
  
  # Determine number of observations in each treatment category
  Z_i <- sparse.model.matrix(~ -1 + factor(Y_i))
  N_t <- apply(Z_i, 2, sum)
  
  # Compute the Average Standardized Absolute Mean (ASAM) difference in the risk factors
  w_mu_t <- t(Z_i / preds) %*% vars / apply(Z_i / preds, 2, sum)
  p_mu <- 1 / N * apply(vars, 2, sum)
  p_sd <- sqrt(diag(t(vars - p_mu) %*% (vars - p_mu)) / (N - 1))
  ASAM <- 1 / T_max * sum(1 / K * apply(t(abs(t(w_mu_t) - p_mu) / p_sd), 1, sum))
  names(ASAM) <- 'ASAM'
  
  # Return the output
  return(c(mnLogLoss(data = data, lev = lev, model = model), ASAM))
}

### Settings for neural network
# General settings
PS_Iter <- 100
seed <- 2019
# Grid search settings
PS_size <- c(1, seq(from = 2, to = 10, by = 2), 25, 50)
PS_decay <- c(0, 0.1, 1, 10, 100, seq(from = 250, to = 1000, by = 50))
PS_maxit <- c(1:20, 30, 40, 50, 75, seq(from = 100, to = 5000, by = 50)) # 100
# Optimal grid
PS_Grid_Opt <- expand.grid(
  size = NA,
  decay = NA,
  maxiter = NA,
  ASAM = NA
)
# Grid 1 - Network-specific parameters
PS_Grid_1 <- expand.grid(
  size = PS_size,
  decay = PS_decay,
  maxiter = PS_Iter,
  ASAM = NA
)
# Grid 2 - Number of iterations
PS_Grid_2 <- expand.grid(
  size = NA,
  decay = NA,
  maxiter = PS_maxit,
  ASAM = NA
)

### Grid search 1
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_NN, selectionFunction = 'best',
                        classProbs = TRUE)
for (i in 1:nrow(PS_Grid_1)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  PS_NN <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                   Undershooting_2 + Premium_New_Base, 
                 data = train_Alt, method = 'nnet', metric = 'logLoss',
                 maximize = FALSE, trControl = cv_ctrl,
                 tuneGrid = data.frame(size = PS_Grid_1$size[i],
                                       decay = PS_Grid_1$decay[i]),
                 trace = FALSE, Hess = FALSE, abstol = 0, reltol = 0,
                 maxit = PS_Grid_1$maxit[i], seed = seed)
  # Store the minimal ASAM obtained
  PS_Grid_1$ASAM[i] <- PS_NN$results$ASAM
}
# Store the optimal parameters from the grid
PS_Grid_Opt$size <- PS_Grid_1$size[which.min(PS_Grid_1$ASAM)]
PS_Grid_Opt$decay <- PS_Grid_1$decay[which.min(PS_Grid_1$ASAM)]
# Adopt these parameters in the next grid searches
PS_Grid_2$size <- PS_Grid_Opt$size
PS_Grid_2$decay <- PS_Grid_Opt$decay
# Determine which iteration in the next grid corresponds to the optimal parameters
Excl_PS_2 <- which(apply(sapply(1:(ncol(PS_Grid_Opt) - 1), function(j) 
  PS_Grid_2[, j] == PS_Grid_1[which.min(PS_Grid_1$ASAM), j]), 1, sum) == (ncol(PS_Grid_Opt) - 1))
# Store these results for that iteration (as they will return the same outcome)
PS_Grid_2$ASAM[Excl_PS_2] <- PS_Grid_1$ASAM[which.min(PS_Grid_1$ASAM)]

### Grid search 2
# Set the cross-validation settings
cv_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_NN, selectionFunction = 'best',
                        classProbs = TRUE)
for (i in c(1:nrow(PS_Grid_2))[-Excl_PS_2]) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Perform cross-validation on the random forest model
  PS_NN <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                   Undershooting_2 + Premium_New_Base, 
                 data = train_Alt, method = 'nnet', metric = 'logLoss',
                 maximize = FALSE, trControl = cv_ctrl,
                 tuneGrid = data.frame(size = PS_Grid_2$size[i],
                                       decay = PS_Grid_2$decay[i]),
                 trace = FALSE, Hess = FALSE, abstol = 0, reltol = 0,
                 maxit = PS_Grid_2$maxit[i], seed = seed)
  # Store the minimal ASAM obtained
  PS_Grid_2$ASAM[i] <- PS_NN$results$ASAM
}
# Store the optimal parameters from the grid
PS_Grid_Opt$maxit <- PS_Grid_2$maxit[which.min(PS_Grid_2$ASAM)]
PS_Grid_Opt$ASAM <- PS_Grid_2$ASAM[which.min(PS_Grid_2$ASAM)]

### Form the cross-validated neural network model
# Set the cross-validation settings
PS_ctrl <- trainControl(method = 'cv', number = 10, search = 'grid',
                        summaryFunction = ASAM_feval_NN, selectionFunction = 'best',
                        classProbs = TRUE)
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Perform cross-validation on the random forest model
PS_NN <- train(form = Offer ~ Policy_Type + Risk_Group + Undershooting_1 + 
                 Undershooting_2 + Premium_New_Base, 
               data = train_Alt, method = 'nnet', metric = 'logLoss',
               maximize = FALSE, trControl = cv_ctrl,
               tuneGrid = data.frame(size = PS_Grid_Opt$size,
                                     decay = PS_Grid_Opt$decay),
               trace = FALSE, Hess = FALSE, abstol = 0, reltol = 0,
               maxit = PS_Grid_Opt$maxit, seed = seed)

### Propensity score estimates
Prop_Score_NN <- predict(PS_NN, newdata = train_Alt, type = 'prob')

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / Prop_Score_NN) %*% Var_XGB / apply(Z / Prop_Score_NN, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_NN <- 1 / K * apply(w_ASAM_t, 2, sum)

##### GAM #####
### Form the Generalized Additive Model
PS_GAM <- gam(formula = list(Offer ~ 1, ~ 1, ~ 1, ~ 1, 1 + 2 + 3 + 4 ~ 
                               s(Undershooting_1, by = Policy_Type, bs = 'cr') + 
                               s(Undershooting_2, by = Policy_Type, bs = 'cr') + 
                               s(Premium_New_Base, by = Policy_Type, bs = 'cr') +
                               s(Undershooting_1, by = Risk_Group, bs = 'cr') + 
                               s(Undershooting_2, by = Risk_Group, bs = 'cr') + 
                               s(Premium_New_Base, by = Risk_Group, bs = 'cr')),
              data = train_Alt, family = multinom(K = T_max - 1), 
              control = gam.control(maxit = 5000, trace = FALSE))

### Propensity score estimates
Prop_Score_GAM <- predict(PS_GAM, newdata = train_Alt, type = 'response')

### Balance after matching
w_mu_t <- as.matrix(t(t(Z / GPS_Score_GAM) %*% Var_XGB / apply(Z / GPS_Score_GAM, 2, sum)))
w_ASAM_t <- abs(w_mu_t - p_mu) / p_sd
Avg_w_ASAM_t_GAM <- 1 / K * apply(w_ASAM_t, 2, sum)

### Store the final results
