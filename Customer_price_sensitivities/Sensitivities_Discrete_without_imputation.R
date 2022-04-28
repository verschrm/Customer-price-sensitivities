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

### Counterfactual response from first imputation only
AMI_Augm$Churn <- as.vector(t(Augm_Response[, , 1]))

##### LASSO #####
### Data for Least Absolute Shrinkage and Selection Operator (LASSO)
# Second order and level interactions
LASSO_df <- data.frame(AMI_Augm$Policy_Type, AMI_Augm$Risk_Group,
                       AMI_Augm$Undershooting_1, AMI_Augm$Undershooting_1 ^ 2,
                       AMI_Augm$Undershooting_2, AMI_Augm$Undershooting_2 ^ 2,
                       AMI_Augm$Premium_New_Base, AMI_Augm$Premium_New_Base ^ 2,
                       AMI_Augm$Competitiveness, AMI_Augm$Competitiveness ^ 2,
                       AMI_Augm$Offer)
LASSO_Cov <- sparse.model.matrix(~ 1 + .^2, data = LASSO_df)
# Exclude third level interactions
Inter_Exc <- c('AMI_Augm.Undershooting_1:AMI_Augm.Undershooting_1.2',
               'AMI_Augm.Undershooting_2:AMI_Augm.Undershooting_2.2',
               'AMI_Augm.Premium_New_Base:AMI_Augm.Premium_New_Base.2',
               'AMI_Augm.Competitiveness:AMI_Augm.Competitiveness.2')
LASSO_Cov <- LASSO_Cov[, which(!(grepl(paste0(Inter_Exc, collapse = '|'), colnames(LASSO_Cov))))]
# Include only interactions with the treatment category and competitiveness
Inter_Incl <- c('AMI_Augm.Offer1', 'AMI_Augm.Offer2', 'AMI_Augm.Offer4', 'AMI_Augm.Offer5',
                'AMI_Augm.Competitiveness', 'AMI_Augm.Competitiveness.2')
LASSO_Cov <- LASSO_Cov[, c(1:18, 18 + which(grepl(paste0(Inter_Incl, collapse = '|'),
                                                  colnames(LASSO_Cov)[19:ncol(LASSO_Cov)])))]

### Force base variables to be used without shrinkage
Pen_Cov <- 1 * (!(colnames(LASSO_Cov) %in%
                    c('(Intercept)', 'AMI_Augm.Premium_New_Base', 'AMI_Augm.Undershooting_1',
                      'AMI_Augm.Undershooting_2', 'AMI_Augm.Risk_Group8ab', 
                      'AMI_Augm.Risk_Group5-7', 'AMI_Augm.Risk_Group1-4', 
                      'AMI_Augm.Policy_TypeEmployee', 'AMI_Augm.Policy_Type2ndCar', 
                      'AMI_Augm.Competitiveness', 'AMI_Augm.Offer1', 
                      'AMI_Augm.Offer2', 'AMI_Augm.Offer4', 'AMI_Augm.Offer5')))

### LASSO parameter selection with cross-validation
# Estimate LASSO model
GLM_LASSO <- cv.glmnet(x = LASSO_Cov, y = AMI_Augm$Churn, family = 'binomial', 
                       type.measure = 'deviance', alpha = 1, nfolds = 10, 
                       lambda = exp(c(seq(from = 5, to = -20, by = -0.01))),
                       penalty.factor = Pen_Cov)
# Retrieve selected parameters
# From: LASSO_Cov[, which(coef(GLM_LASSO, s = 'lambda.1se')[-1] != 0)]
MI_formula <- as.formula('Churn ~ 1 + Policy_Type + Risk_Group + Undershooting_1 + Undershooting_2 +
                                       Premium_New_Base + Competitiveness + Offer +
                         I(Undershooting_1 ^ 2) + I(Premium_New_Base ^ 2) +
                         I(1 * (Policy_Type == "Employee")) : Competitiveness +
                         I(1 * (Policy_Type == "Employee")) : I(Competitiveness ^ 2) +
                         I(1 * (Policy_Type == "Employee")) : I(1 * (Offer == "1")) +
                         I(1 * (Policy_Type == "2ndCar")) : I(1 * (Offer == "1")) +
                         I(1 * (Policy_Type == "2ndCar")) : I(1 * (Offer == "2")) +
                         I(1 * (Policy_Type == "2ndCar")) : I(1 * (Offer == "4")) +
                         I(1 * (Policy_Type == "2ndCar")) : I(1 * (Offer == "5")) +
                         I(1 * (Risk_Group == "8ab")) : Competitiveness + 
                         I(1 * (Risk_Group == "8ab")) : I(1 * (Offer == "1")) +
                         I(1 * (Risk_Group == "8ab")) : I(1 * (Offer == "5")) +
                         I(1 * (Risk_Group == "5-7")) : Competitiveness + 
                         I(1 * (Risk_Group == "5-7")) : I(1 * (Offer == "1")) +
                         I(1 * (Risk_Group == "5-7")) : I(1 * (Offer == "2")) +
                         I(1 * (Risk_Group == "5-7")) : I(1 * (Offer == "4")) +
                         I(1 * (Risk_Group == "5-7")) : I(1 * (Offer == "5")) +
                         I(1 * (Risk_Group == "1-4")) : I(1 * (Offer == "5")) +
                         Undershooting_1 : I(1 * (Offer == "2")) +
                         Undershooting_1 : I(1 * (Offer == "5")) +
                         I(Undershooting_1 ^ 2) : I(1 * (Offer == "4")) +
                         Undershooting_2 : Competitiveness + Undershooting_2 : I(1 * (Offer == "1")) + 
                         Undershooting_2 : I(1 * (Offer == "2")) + 
                         Undershooting_2 : I(1 * (Offer == "4")) + 
                         I(Undershooting_2 ^ 2) : I(1 * (Offer == "1")) +
                         I(Undershooting_2 ^ 2) : I(1 * (Offer == "2")) +
                         I(Undershooting_2 ^ 2) : I(1 * (Offer == "5")) +
                         Premium_New_Base : Competitiveness +
                         Premium_New_Base : I(1 * (Offer == "1")) + 
                         Premium_New_Base : I(1 * (Offer == "4")) + 
                         Premium_New_Base : I(1 * (Offer == "5")) + 
                         I(Premium_New_Base ^ 2) : I(1 * (Offer == "1")) + 
                         Competitiveness : I(1 * (Offer == "1")) + 
                         Competitiveness : I(1 * (Offer == "4")) + 
                         I(Competitiveness ^ 2) : I(1 * (Offer == "1")) + 
                         I(Competitiveness ^ 2) : I(1 * (Offer == "4")) + 
                         I(Competitiveness ^ 2): I(1 * (Offer == "5"))')

##### Causal inference #####
### Discrete global response model
Logit_GLM <- glm(MI_formula, data = AMI, family = binomial(link = 'logit'))

### Churn predictions
Churn_Pred <- predict(Logit_GLM, AMI_Augm, type = 'response')
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
  Churn <- predict(Churn_Model, Data_Row, type = 'response')
  
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
    sum(Churn_Pred(Churn_Model = Logit_GLM,
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
    sum(Churn_Pred(Churn_Model = Logit_GLM, Offers = as.character(t), 
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
    sum(Churn_Pred(Churn_Model = Logit_GLM, Offers = AMI_r$Offer,
                   Comp = C_grid[c], Data_Row = AMI_r)))
}

### Store the final results
