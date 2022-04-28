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
Packages <- c('dplyr', 'fastDummies', 'Matrix', 'lme4', 'gbm', 'xgboost', 'Rsymphony')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load output from global response model procedure
# Load results from Sensitivities_Discrete_XGBoost.R

##### LP optimization #####
### Churn prediction function for a single period
Churn_Pred_i_1 <- function(Churn_Model, Cat_Offer, Offers, Data_i_1) {
  # Offers to consider
  Data_i_1$Offer[] <- Cat_Offer
  
  # Determine churn predictions
  Var_1_XGB <- as.matrix(dummy_cols(Data_i_1[, c('Churn', 'Competitiveness', 'Offer',
                                                 'Policy_Type', 'Risk_Group', 'Undershooting_1', 
                                                 'Undershooting_2', 'Premium_New_Base')], 
                                    remove_most_frequent_dummy = FALSE)[, -c(1, 3, 4, 5)])
  Var_Row_XGB <- Matrix(Var_Row_XGB, sparse = TRUE)
  Churn_i_1 <- predict(Churn_Model, newdata = as.matrix(Var_1_XGB), type = 'response')
  
  # Return the output
  return(Churn_i_1)
}

### LP optimization function for a single period
LP_Rsymphony_optim_i_1 <- function(Churn_i_1, P_Old_i_1, C_i_1, RHS_con, Sign_con, Sparse_con, 
                                   alpha, Offers) {
  # Preliminary variables
  T_max <- length(Offers)
  N_1 <- ncol(Churn_i_1)
  
  # Objective coefficients of LP, where we can eliminate one decision variable per renewal
  # since the sum of the decision variables needs to equal one for each renewal
  LP_obj <- t(1 - Churn_i_1) * (P_Old_i_1 %o% (1 + Offers) - C_i_1)
  LP_obj <- as.vector(t(LP_obj[, -T_max] - LP_obj[, T_max]))
  
  # Constraint coefficients and signs of LP
  LP_RHS <- Churn_i_1[T_max, ]
  LP_con_sparse <- rbind(Sparse_con[1:N_1, ], as.vector(t(t(Churn_i_1[-T_max, ]) - LP_RHS)) / N_1)
  LP_RHS <- c(head(RHS_con, -1), alpha - sum(LP_RHS) / N_1)
  LP_sign <- Sign_con
  
  # Optimization of LP
  LP_opt <- Rsymphony_solve_LP(max = TRUE, obj = LP_obj, dir = LP_sign, rhs = LP_RHS, 
                               mat = LP_con_sparse, time_limit = 3600,
                               types = rep('B', length(LP_obj)))
  
  # Optimal offers
  LP_Ind <- sapply(1:N_1, function(i) min(T_max, which(LP_opt$solution[
    (1 + (T_max - 1) * (i - 1)):((T_max - 1) * i)] == 1)))
  
  # Return the output
  return(rbind(LP_opt$status, LP_Ind))
}

### Preliminary variables
N_Policies <- nrow(AMI)
Q_means <- sapply(1:T_max, function(t) (Q_bins[t] + Q_bins[t + 1]) / 2)
Offer_Grid <- as.numeric(Q_means)
P_Old <- AMI$Premium_Old
C_1 <- AMI$Costs

### Churn predictions for a single period
Churn_i <- t(sapply(1:length(Offer_Grid), function(t_1)
  Churn_Pred_i_1(Churn_Model = Churn_XGB, Cat_Offer = as.character(t_1), Data_i_1 = AMI)))
Profit_i <- (1 - Churn_i) * t(P_Old %o% (1 + Offer_Grid) - C_1)
alpha_Range <- apply(sapply(1:N_Policies, function(i) 
  c(min(Churn_i[, i]), Churn_i[which.max(Profit_i[, i]), i])), 1, mean)

### LP setup
# Seed for the random number generator
seed <- 20200101
# Sparse matrix of nrow constraints and ncol variables
LP_con_sparse <- Matrix(0, nrow = N_Policies + 1, ncol = N_Policies * (length(Offer_Grid) - 1), 
                        byrow = TRUE, sparse = TRUE)
LP_con_sparse[cbind(rep(1:N_Policies, each = length(Offer_Grid) - 1), 
                    1:(N_Policies * (length(Offer_Grid) - 1)))] <- 1
LP_con_sparse[N_Policies + 1, ] <- 1
# Constraints' signs
LP_sign <- c(rep('<=', N_Policies), '<=')
# Maximal admissible overall portfolio churn rate
a_Step <- 0.05 / 100
alpha_Grid <- c(alpha_Range[1], seq(from = ceiling(alpha_Range[1] / a_Step) * a_Step, 
                                    to = floor(alpha_Range[2] / a_Step) * a_Step, by = a_Step), 
                alpha_Range[2], 0.50)
# Constraints' RHSs
LP_RHS <- c(rep(1, N_Policies), NA)

### Pre-allocation of output
Results_LP <- array(NA, dim = c(length(alpha_Grid), 5, N_Policies))

### Boundary solutions
Offers_min <- sapply(1:N_Policies, function(i) which.min(Churn_i[, i]))
Offers_max <- sapply(1:N_Policies, function(i) which.max(Profit_i[, i]))

### Loop over the alpha grid of admissible overall portfolio churn rates
for (a in 1:length(alpha_Grid)) {
  # For first alpha use boundary solution
  if (a == 1) {
    Results_LP[a, 1:2, ] <- rbind(0, Offers_min)
    # For final two alpha's use boundary solution
  } else if (a >= (length(alpha_Grid) - 1)) {
    Results_LP[a, 1:2, ] <- rbind(0, Offers_max)
  } else {
    # Ensure reproducibility by setting the seed for the random number generator
    set.seed(seed)
    # Return status of solution and offer indices
    Results_LP[a, 1:2, ] <- t(matrix(LP_Rsymphony_optim_i_1(Churn_i_1 = Churn_i, P_Old_i_1 = P_Old, 
                                                            C_i_1 = C_1, RHS_con = LP_RHS, 
                                                            Sign_con = LP_sign, 
                                                            Sparse_con = LP_con_sparse, 
                                                            alpha = alpha_Grid[a], 
                                                            Offers = Offer_Grid), 
                                     ncol = 2, byrow = TRUE))
  }
  # In case no optimal solution has been reached, use solution to previous alpha
  if (Results_LP[a, 1, 1] > 0) {
    Results_LP[a, 3, ] <- Churn_i[cbind(Results_LP[a - 1, 2, ], 1:N_Policies)]
    Results_LP[a, 4, ] <- Offer_Grid[Results_LP[a - 1, 2, ]]
    # Otherwise retrieve corresponding churn estimates and renewal offers
  } else {
    Results_LP[a, 3, ] <- Churn_i[cbind(Results_LP[a, 2, ], 1:N_Policies)]
    Results_LP[a, 4, ] <- Offer_Grid[Results_LP[a, 2, ]]
  }
  # Calculate expected optimal profit
  Results_LP[a, 5, ] <- (1 - Results_LP[a, 3, ]) * (P_Old * (1 + Results_LP[a, 4, ]) - C_1)
}

##### Efficient frontier #####
### Actual profit and churn
Act <- cbind((1 - AMI$Churn) %*% (AMI$Premium_New - C_1), mean(AMI$Churn))

### Optimal profit and churn
Opt <- cbind(as.numeric(alpha_Grid),
             apply(Results_LP[, 5, ], 1, sum),
             apply(Results_LP[, 3, ], 1, mean))

### Relative profit and churn
Rel <- cbind(Opt[, 1], Opt[, 2] / Act[, 1] - 1, Opt[, 3], (Opt[, 2] - Act[, 1]) / N_Policies)

### Store the final results
