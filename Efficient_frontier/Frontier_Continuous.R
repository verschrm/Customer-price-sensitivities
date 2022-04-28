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
Packages <- c('dplyr', 'fastDummies', 'Matrix', 'lme4', 'gbm', 'xgboost', 'Rsymphony', 'nloptr')
invisible(lapply(Packages, library, character.only = TRUE))

rm(Packages)

##### Load data #####
### Load output from conditional dose-response function procedure
# Load results from Sensitivities_Continuous.R

##### (N)LP settings #####
### Churn prediction function for a single period
Churn_Pred_i_1 <- function(Churn_Model, Offers, GPS_Model, sd_GPS, sd_t, Data_i_1) {
  # Preliminary variables
  T_t <- length(Offers)
  N_1 <- nrow(Data_i_1)
  
  # Generalized propensity score input
  Var_i_XGB <- as.matrix(dummy_cols(Data_i_1[, c('Policy_Type', 'Risk_Group', 'Undershooting_1',
                                                 'Undershooting_2', 'Premium_New_Base')], 
                                    remove_most_frequent_dummy = FALSE)[, -c(1, 2)])
  Var_i_XGB <- Matrix(Var_i_XGB, sparse = TRUE)
  Mean <- predict(GPS_Model, newdata = as.matrix(Var_i_XGB), type = 'response')
  Normalizer <- apply(sapply(1:(length(Q_bins) - 1), function(t) 
    pnorm((Q_bins[t + 1] - Mean) / sd_GPS[t], mean = 0, sd = 1) - 
      pnorm((Q_bins[t] - Mean) / sd_GPS[t], mean = 0, sd = 1)), 1, sum)
  
  # Offers to consider
  Data_i_1 <- Data_i_1[rep(1:N_1, each = T_t), ]
  Data_i_1$Offer <- rep(Offers, N_1)
  
  # Estimated GPS
  Data_i_1$GPS <- as.vector(t(sapply(1:T_t, function(t)
    dnorm((Offers[t] - Mean) / sd_t[t], mean = 0, sd = 1) / Normalizer / sd_t[t])))
  
  # Determine churn predictions
  Var_i_XGB <- as.matrix(dummy_cols(
    Data_i_1[, c('Competitiveness', 'Offer', 'GPS', 'Policy_Type', 'Risk_Group',
                 'Undershooting_1', 'Undershooting_2', 'Premium_New_Base')],
    remove_most_frequent_dummy = FALSE)[, -c(4, 5)])
  Var_i_XGB <- Matrix(Var_i_XGB, sparse = TRUE)
  Churn_i_1 <- predict(Churn_Model, newdata = as.matrix(Var_i_XGB), type = 'response')
  Churn_i_1 <- matrix(Churn_i_1, nrow = N_1, ncol = T_t, byrow = TRUE)
  
  # Return the output
  return(Churn_i_1)
}

### NLP constrained objective function for a single period
NLP_obj_con_i_1 <- function(Offers, sd_GPS, Data_i_1, Mean, Normalizer, Churn_Model, 
                              P_Old_i_1, C_i_1, eps, alpha) {
  # Offers to consider
  Data_i_1$Offer <- Offers
  sd_t <- sd_GPS[as.numeric(cut(Offers, breaks = Q_bins, include.lowest = TRUE))]
  
  # Estimated GPS
  Data_i_1$GPS <- dnorm((Offers - Mean) / sd_t, mean = 0, sd = 1) / Normalizer / sd_t
  
  # Estimated churn rate
  Var_i_XGB <- as.matrix(dummy_cols(
    Data_i_1[, c('Competitiveness', 'Offer', 'GPS', 'Policy_Type', 'Risk_Group',
                 'Undershooting_1', 'Undershooting_2', 'Premium_New_Base')],
    remove_most_frequent_dummy = FALSE)[, -c(4, 5)])
  Var_i_XGB <- Matrix(Var_i_XGB, sparse = TRUE)
  Churn_i_1 <- predict(Churn_Model, newdata = as.matrix(Var_i_XGB), type = 'response')
  
  # Objective value of NLP
  NLP_obj <- -(1 - Churn_i_1) %*% (P_Old_i_1 * (1 + Offers) - C_i_1)
  
  # Check if the inequality constraint is met
  NLP_con <- sum(Churn_i_1) / nrow(Data_i_1) - alpha
  NLP_obj_con <- (NLP_con <= eps) * NLP_obj
  
  # Return the output
  return(NLP_obj_con)
}

### NLP churn estimates for a single period
NLP_Churn_i_1 <- function(Offers, sd_GPS, Data_i_1, Mean, Normalizer, Churn_Model) {
  # Offers to consider
  Data_i_1$Offer <- Offers
  sd_t <- sd_GPS[as.numeric(cut(Offers, breaks = Q_bins, include.lowest = TRUE))]
  
  # Estimated GPS
  Data_i_1$GPS <- dnorm((Offers - Mean) / sd_t, mean = 0, sd = 1) / Normalizer / sd_t
  
  # Estimated churn rate
  Var_i_XGB <- as.matrix(dummy_cols(
    Data_i_1[, c('Competitiveness', 'Offer', 'GPS', 'Policy_Type', 'Risk_Group',
                 'Undershooting_1', 'Undershooting_2', 'Premium_New_Base')],
    remove_most_frequent_dummy = FALSE)[, -c(4, 5)])
  Var_i_XGB <- Matrix(Var_i_XGB, sparse = TRUE)
  Churn_i_1 <- predict(Churn_Model, newdata = as.matrix(Var_i_XGB), type = 'response')
  
  # Return the output
  return(Churn_i_1)
}

### LP local optimization function for a single period
LP_L_Rsymphony_optim_i_1 <- function(Data_i_1, sd_GPS, Mean, Normalizer, Churn_Model,
                                     RHS_con, Sign_con, Sparse_con, alpha, 
                                     Offers_Init, Num, Scale, LB, UB) {
  # Preliminary variables
  N_1 <- length(Offers_Init)
  Scale_Num <- Scale / (Num - 1)
  
  # Determine the local search domain
  Scale_L_temp <- ifelse((Offers_Init - Scale / 2) < LB, Offers_Init - LB, Scale / 2)
  Scale_U_temp <- ifelse((Offers_Init + Scale / 2) > UB, UB - Offers_Init, Scale / 2)
  Scale_L <- Scale_L_temp + (Scale / 2 - Scale_U_temp)
  Scale_U <- Scale_U_temp + (Scale / 2 - Scale_L_temp)
  Offers <- mapply(function(t, l, u) 
    c(t - l, seq(from = t - ceiling(l / Scale_Num - 1) * Scale_Num, 
                 to = t + floor(u / Scale_Num - 1) * Scale_Num,
                 length.out = Num - 2), t + u), Offers_Init, Scale_L, Scale_U)
  
  # Estimated churn rate
  Churn_i_1 <- t(sapply(1:Num, function(t) 
    NLP_Churn_i_1(Offers = Offers[t, ], sd_GPS = sd_GPS, Data_i_1 = Data_i_1, 
                  Mean = Mean, Normalizer = Normalizer, Churn_Model = Churn_Model)))
  
  # Objective coefficients of LP, where we can eliminate one decision variable per renewal
  # since the sum of the decision variables needs to equal one for each renewal
  LP_Profit <- t(1 - Churn_i_1) * (Data_i_1$Premium_Old * t(1 + Offers) - Data_i_1$Costs)
  LP_obj <- -as.vector(t(LP_Profit[, -Num] - LP_Profit[, Num]))
  
  # Constraint coefficients and signs of LP
  LP_RHS <- Churn_i_1[Num, ]
  LP_con_sparse <- rbind(Sparse_con[1:N_1, ], as.vector(t(t(Churn_i_1[-Num, ]) - LP_RHS)) / N_1)
  LP_RHS <- c(head(RHS_con, -1), alpha - sum(LP_RHS) / N_1)
  LP_sign <- Sign_con
  
  # Optimization of LP
  LP_opt <- Rsymphony_solve_LP(max = FALSE, obj = LP_obj, dir = LP_sign, rhs = LP_RHS, 
                               mat = LP_con_sparse, types = rep('B', length(LP_obj)),
                               time_limit = 600,
                               verbosity = -2)
  
  # Optimal offers
  LP_Ind <- sapply(1:N_1, function(i) min(Num, which(LP_opt$solution[
    (1 + (Num - 1) * (i - 1)):((Num - 1) * i)] == 1)))
  LP_Offers <- Offers[cbind(LP_Ind, 1:N_1)]
  
  # Optimal churn and profit
  LP_Churn <- Churn_i_1[cbind(LP_Ind, 1:N_1)]
  LP_Profit <- LP_Profit[cbind(1:N_1, LP_Ind)]
  
  # Return the output
  return(rbind(LP_opt$status, LP_Offers, LP_Churn, LP_Profit))
}

### Preliminary variables
N_Policies <- nrow(AMI)
P_Old <- AMI$Premium_Old
C_1 <- AMI$Costs

### Generalized propensity score input for a single period
Var_i_XGB_1 <- as.matrix(dummy_cols(AMI[, c('Policy_Type', 'Risk_Group', 'Undershooting_1',
                                            'Undershooting_2', 'Premium_New_Base')], 
                                    remove_most_frequent_dummy = FALSE)[, -c(1, 2)])
Var_i_XGB_1 <- Matrix(Var_i_XGB_1, sparse = TRUE)
Mean_1 <- predict(GPS_XGB, newdata = as.matrix(Var_i_XGB_1), type = 'response')
Normalizer_1 <- apply(sapply(1:(length(Q_bins) - 1), function(t) 
  pnorm((Q_bins[t + 1] - Mean_1) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Mean_1) / sd_GPS_XGB[t], mean = 0, sd = 1)), 1, sum)

##### Initial local search #####
### Initial local search options
Step_Shift <- 1/10 / 100
Tries <- ceiling((max(T_Offers) - min(T_Offers) - Step_Shift) / Step_Shift)
Step_Num <- 1000
ILS_Tries <- array(NA, dim = c(2 * Tries + 1, 2, 2))
ILS_Offers <- array(NA, dim = c(2 * Tries + 1, N_Policies, 2))

### Initial local search to find minimal churn and maximal profit offers
for (try in (-Tries):Tries) {
  # Initial local search offers
  Offers_try <- seq(from = min(T_Offers) + max(try, 0) * Step_Shift, 
                    to = max(T_Offers) + min(try, 0) * Step_Shift, length.out = Step_Num)
  sd_try <- sd_GPS_XGB[as.numeric(cut(Offers_try, breaks = Q_bins, include.lowest = TRUE))]
  
  # Estimated churn rate
  Churn_try <- t(Churn_Pred_i_1(Churn_Model = Churn_XGB, Offers = Offers_try, 
                                GPS_Model = GPS_XGB, sd_GPS = sd_GPS_XGB, sd_t = sd_try, 
                                Data_i_1 = AMI))
  
  # Predicted profit
  Profit_try <- (1 - Churn_try) * t(P_Old %o% (1 + Offers_try) - C_1)
  
  # Minimal and maximal indices
  Ind_try <- sapply(1:N_Policies, function(i) 
    c(which.min(Churn_try[, i]), which.max(Profit_try[, i])))
  
  # Corresponding portfolio churn rate and profit
  ILS_Tries[try + Tries + 1, , ] <- cbind(
    rbind(sum(Churn_try[cbind(Ind_try[1, ], 1:N_Policies)]) / N_Policies,
          sum(Profit_try[cbind(Ind_try[1, ], 1:N_Policies)])),
    rbind(sum(Churn_try[cbind(Ind_try[2, ], 1:N_Policies)]) / N_Policies,
          sum(Profit_try[cbind(Ind_try[2, ], 1:N_Policies)])))
  
  # Corresponding offers
  ILS_Offers[try + Tries + 1, , ] <- cbind(Offers_try[Ind_try[1, ]], Offers_try[Ind_try[2, ]])
}

### Tabulate results
# Actual expected churn and profit
Churn_Act <- NLP_Churn_i_1(Offers = AMI$Offer, sd_GPS = sd_GPS_XGB, Data_i_1 = AMI, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
Profit_Act <- as.numeric((1 - Churn_Act) %*% (P_Old * (1 + AMI$Offer) - C_1))
# Table format of initial local search results
ILS_df <- data.frame(cbind((-Tries):Tries, ILS_Tries[, 1, ], ILS_Tries[, 2, ] / Profit_Act - 1))
colnames(ILS_df) <- c('Iteration', 'Churn_min', 'Churn_opt', 'Profit_min', 'Profit_opt')
# Restructure table format
ILS_Iter_df <- data.frame(rbind(cbind(ILS_df$Iteration, ILS_df$Churn_min, ILS_df$Profit_min, 
                                     (ILS_df$Profit_min + 1) * Profit_Act),
                               cbind(ILS_df$Iteration, ILS_df$Churn_opt, ILS_df$Profit_opt, 
                                     (ILS_df$Profit_opt + 1) * Profit_Act)))
colnames(ILS_Iter_df) <- c('Iteration', 'Churn_Avg', 'Profit_Tot_Rel', 'Profit_Tot_Abs')
ILS_Init_df <- as.matrix(rbind(ILS_Offers[, , 1], ILS_Offers[, , 2]))

##### (N)LP optimization #####
### Maximal admissible overall portfolio churn rate
a_Step <- 0.1 / 100
a_Focus <- c(ceiling(min(ILS_Iter_df$Churn_Avg) / a_Step) * a_Step, 0.50)
alpha_Grid <- seq(from = ceiling(a_Focus[1] / a_Step) * a_Step,
                  to = floor(a_Focus[2] / a_Step) * a_Step, by = a_Step)

### (N)LP setup
# Seed for the random number generator
seed <- 20200101
# LP local search domain properties
Num_L <- 3
Scale_1 <- 2 * 2 / 100
Scale_2 <- 2 * 1 / 100
Scale_3 <- 2 * 0.5 / 100
Scale_4 <- 2 * 0.25 / 100
Scale_5 <- 2 * 0.125 / 100
# LP sparse matrix of nrow constraints and ncol variables
LP_con_sparse_L <- Matrix(0, nrow = N_Policies, ncol = N_Policies * (Num_L - 1), 
                          byrow = TRUE, sparse = TRUE)
LP_con_sparse_L[cbind(rep(1:N_Policies, each = Num_L - 1), 1:(N_Policies * (Num_L - 1)))] <- 1
LP_con_sparse_L <- rbind(LP_con_sparse_L, 1)
# LP constraints' signs
LP_sign <- c(rep('<=', N_Policies), '<=')
# LP constraints' RHSs
LP_RHS <- c(rep(1, N_Policies), NA)
# (N)LP constraint tolerance
LP_eps <- 0
NLP_eps <- 1e-6
# NLP optimization algorithm
NLP_Opts <- list('algorithm' = 'NLOPT_LN_SBPLX', 'ftol_rel' = 0, 'ftol_abs' = 0,
                 'xtol_rel' = 0, 'xtol_abs' = rep(0, N_Policies),
                 'maxeval' = 1000, 'print_level' = 0)
# NLP lower and upper bounds of decision variables
NLP_LB <- rep(as.numeric(Q_bins[1]), N_Policies)
NLP_UB <- rep(as.numeric(tail(Q_bins, 1)), N_Policies)

### Pre-allocation of output
Results_LP_1 <- array(NA, dim = c(length(alpha_Grid), 4, N_Policies))
Results_LP_2 <- array(NA, dim = c(length(alpha_Grid), 4, N_Policies))
Results_LP_3 <- array(NA, dim = c(length(alpha_Grid), 4, N_Policies))
Results_LP_4 <- array(NA, dim = c(length(alpha_Grid), 4, N_Policies))
Results_LP_5 <- array(NA, dim = c(length(alpha_Grid), 4, N_Policies))
Results_NLP <- array(NA, dim = c(length(alpha_Grid), 4, N_Policies))

### Initial values
Ind_0 <- 0
NLP_Init_a <- rep(NA, N_Policies)
Init_0 <- rep(0, N_Policies)
Profit_a <- 0
NLP_Profit_a <- 0

### Loop over the alpha grid of admissible overall portfolio churn rates
for (a in 1:length(alpha_Grid)) {
  # Select closest, optimal initial solution from initial local search
  Iter_a <- which(ILS_Iter_df$Churn_Avg <= (alpha_Grid[a] + NLP_eps))
  Ind_a <- Iter_a[which.max(ILS_Iter_df$Profit_Tot_Abs[Iter_a])]
  # Check if this leads to a different starting point
  if (Ind_a != Ind_0) {
    Init_a <- ILS_Init_df[Ind_a, ]
    Profit_a <- -ILS_Iter_df$Profit_Tot_Abs[Ind_a]
  }
  # Check if previous optimal solution leads to a more profitable starting point
  if (NLP_Profit_a < Profit_a) {
    Init_a <- NLP_Init_a
    Profit_a <- NLP_Profit_a
  }
  
  ### Try to improve initial solution by LP local optimization
  # Only initiate LP optimization if starting point is different from previous starting point
  if (!identical(Init_a, Init_0)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Return status of solution and offer indices for scale 1
  Results_LP_1[a, 1:4, ] <- LP_L_Rsymphony_optim_i_1(Data_i_1 = AMI, sd_GPS = sd_GPS_XGB, 
                                                     Mean = Mean_1, Normalizer = Normalizer_1, 
                                                     Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                                                     Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                                                     alpha = alpha_Grid[a] + LP_eps, 
                                                     Offers_Init = Init_a, Num = Num_L, 
                                                     Scale = Scale_1, LB = as.numeric(Q_bins[1]), 
                                                     UB = as.numeric(tail(Q_bins, 1)))
  # Adjust starting points if solution is feasible and better
  if (((-sum(Results_LP_1[a, 4, ])) < Profit_a) & 
      ((floor(mean(Results_LP_1[a, 3, ]) / NLP_eps) * NLP_eps) <= (alpha_Grid[a] + NLP_eps))) {
    Init_a <- Results_LP_1[a, 2, ]
    Profit_a <- -sum(Results_LP_1[a, 4, ])
  }
  
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Return status of solution and offer indices for scale 2
  Results_LP_2[a, 1:4, ] <- LP_L_Rsymphony_optim_i_1(Data_i_1 = AMI, sd_GPS = sd_GPS_XGB, 
                                                     Mean = Mean_1, Normalizer = Normalizer_1, 
                                                     Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                                                     Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                                                     alpha = alpha_Grid[a] + LP_eps, 
                                                     Offers_Init = Init_a, Num = Num_L, 
                                                     Scale = Scale_2, LB = as.numeric(Q_bins[1]), 
                                                     UB = as.numeric(tail(Q_bins, 1)))
  # Adjust starting points if solution is feasible and better
  if (((-sum(Results_LP_2[a, 4, ])) < Profit_a) & 
      ((floor(mean(Results_LP_2[a, 3, ]) / NLP_eps) * NLP_eps) <= (alpha_Grid[a] + NLP_eps))) {
    Init_a <- Results_LP_2[a, 2, ]
    Profit_a <- -sum(Results_LP_2[a, 4, ])
  }
  
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Return status of solution and offer indices for scale 3
  Results_LP_3[a, 1:4, ] <- LP_L_Rsymphony_optim_i_1(Data_i_1 = AMI, sd_GPS = sd_GPS_XGB, 
                                                     Mean = Mean_1, Normalizer = Normalizer_1, 
                                                     Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                                                     Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                                                     alpha = alpha_Grid[a] + LP_eps, 
                                                     Offers_Init = Init_a, Num = Num_L, 
                                                     Scale = Scale_3, LB = as.numeric(Q_bins[1]), 
                                                     UB = as.numeric(tail(Q_bins, 1)))
  # Adjust starting points if solution is feasible and better
  if (((-sum(Results_LP_3[a, 4, ])) < Profit_a) & 
      ((floor(mean(Results_LP_3[a, 3, ]) / NLP_eps) * NLP_eps) <= (alpha_Grid[a] + NLP_eps))) {
    Init_a <- Results_LP_3[a, 2, ]
    Profit_a <- -sum(Results_LP_3[a, 4, ])
  }
  
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Return status of solution and offer indices for scale 4
  Results_LP_4[a, 1:4, ] <- LP_L_Rsymphony_optim_i_1(Data_i_1 = AMI, sd_GPS = sd_GPS_XGB, 
                                                     Mean = Mean_1, Normalizer = Normalizer_1, 
                                                     Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                                                     Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                                                     alpha = alpha_Grid[a] + LP_eps, 
                                                     Offers_Init = Init_a, Num = Num_L, 
                                                     Scale = Scale_4, LB = as.numeric(Q_bins[1]), 
                                                     UB = as.numeric(tail(Q_bins, 1)))
  # Adjust starting points if solution is feasible and better
  if (((-sum(Results_LP_4[a, 4, ])) < Profit_a) & 
      ((floor(mean(Results_LP_4[a, 3, ]) / NLP_eps) * NLP_eps) <= (alpha_Grid[a] + NLP_eps))) {
    Init_a <- Results_LP_4[a, 2, ]
    Profit_a <- -sum(Results_LP_4[a, 4, ])
  }
  
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Return status of solution and offer indices for scale 5
  Results_LP_5[a, 1:4, ] <- LP_L_Rsymphony_optim_i_1(Data_i_1 = AMI, sd_GPS = sd_GPS_XGB, 
                                                     Mean = Mean_1, Normalizer = Normalizer_1, 
                                                     Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                                                     Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                                                     alpha = alpha_Grid[a] + LP_eps, 
                                                     Offers_Init = Init_a, Num = Num_L, 
                                                     Scale = Scale_5, LB = as.numeric(Q_bins[1]), 
                                                     UB = as.numeric(tail(Q_bins, 1)))
  
  # Adjust starting points if solution is feasible and better
  if (((-sum(Results_LP_5[a, 4, ])) < Profit_a) & 
      ((floor(mean(Results_LP_5[a, 3, ]) / NLP_eps) * NLP_eps) <= (alpha_Grid[a] + NLP_eps))) {
    NLP_Init_a <- Results_LP_5[a, 2, ]
    NLP_Profit_a <- -sum(Results_LP_5[a, 4, ])
  } else {
    NLP_Init_a <- Init_a
    NLP_Profit_a <- Profit_a
  }
  } else {
    # If starting point did not change, use previous solution
    Results_LP_1[a, 1:4, ] <- Results_LP_1[a - 1, 1:4, ]
    Results_LP_2[a, 1:4, ] <- Results_LP_2[a - 1, 1:4, ]
    Results_LP_3[a, 1:4, ] <- Results_LP_3[a - 1, 1:4, ]
    Results_LP_4[a, 1:4, ] <- Results_LP_4[a - 1, 1:4, ]
    Results_LP_5[a, 1:4, ] <- Results_LP_5[a - 1, 1:4, ]
  }
  
  ### Try to improve initial/LP solution by NLP optimization
  # Specify specific objective function with inequality constraint
  NLP_obj_con_a <- function(x) {
    return(NLP_obj_con_i_1(Offers = x, sd_GPS = sd_GPS_XGB, Data_i_1 = AMI, Mean = Mean_1, 
                           Normalizer = Normalizer_1, Churn_Model = Churn_XGB, 
                           P_Old_i_1 = P_Old, C_i_1 = C_1, eps = NLP_eps, 
                           alpha = alpha_Grid[a]))
  }
  # Only initiate NLP optimization if starting point is different from previous starting point
  if (!identical(NLP_Init_a, Init_0)) {
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # Retrieve NLP solution
  NLP_opt <- nloptr(x0 = NLP_Init_a, eval_f = NLP_obj_con_a,
                    lb = NLP_LB, ub = NLP_UB, opts = NLP_Opts)
  # Return status of solution
  Results_NLP[a, 1, ] <- NLP_opt$status
  # Retrieve corresponding churn estimates and renewal offers
  Results_NLP[a, 2, ] <- NLP_opt$solution
  Results_NLP[a, 3, ] <- NLP_Churn_i_1(Offers = Results_NLP[a, 2, ], sd_GPS = sd_GPS_XGB, 
                                       Data_i_1 = AMI, Mean = Mean_1, 
                                       Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
  # Calculate expected optimal profit
  Results_NLP[a, 4, ] <- (1 - Results_NLP[a, 3, ]) * (P_Old * (1 + Results_NLP[a, 2, ]) - C_1)
  # Adjust solution if starting points are better or solution is infeasible
  if ((NLP_Profit_a < (-sum(Results_NLP[a, 4, ]))) |
      ((alpha_Grid[a] + NLP_eps) < (floor(mean(Results_NLP[a, 3, ]) / NLP_eps) * NLP_eps))) {
    # Retrieve corresponding churn estimates and renewal offers
    Results_NLP[a, 2, ] <- NLP_Init_a
    Results_NLP[a, 3, ] <- NLP_Churn_i_1(Offers = Results_NLP[a, 2, ], sd_GPS = sd_GPS_XGB, 
                                         Data_i_1 = AMI, Mean = Mean_1, 
                                         Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
    # Calculate expected optimal profit
    Results_NLP[a, 4, ] <- (1 - Results_NLP[a, 3, ]) * (P_Old * (1 + Results_NLP[a, 2, ]) - C_1)
  }
  } else {
    # If starting point did not change, use previous solution
    Results_NLP[a, 1:4, ] <- Results_NLP[a - 1, 1:4, ]
  }
  
  # Use current solution as initial values for next optimization
  NLP_Init_a <- Results_NLP[a, 2, ]
  NLP_Profit_a <- -sum(Results_NLP[a, 4, ])
  # Allocate previous optimal solution
  Ind_0 <- Ind_a
  if (a > 1) {
    Init_0 <- Results_NLP[a - 1, 2, ]
  }
}

##### Efficient frontier #####
### Actual profit and churn
Act <- cbind((1 - AMI$Churn) %*% (AMI$Premium_New - C_1), mean(AMI$Churn))

### Optimal profit and churn
Opt <- cbind(as.numeric(alpha_Grid),
             apply(Results_NLP[, 4, ], 1, sum),
             apply(Results_NLP[, 3, ], 1, mean))

### Relative profit and churn
Rel <- cbind(Opt[, 1], Opt[, 2] / Act[, 1] - 1, Opt[, 3], (Opt[, 2] - Act[, 1]) / N_Policies)

### Store the final results
