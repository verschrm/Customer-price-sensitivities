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
### Load output from global response model and efficient frontier procedure
# Load results from Frontier_Discrete_without_imputation.R

##### (N)LP settings #####
### Churn prediction function for a single period without feedback adjustments
Churn_Pred_i_1 <- function(Churn_Model, Cat_Offer, Data_i_1) {
  # Offers to consider
  Data_i_1$Offer[] <- Cat_Offer
  
  # Determine churn predictions
  Churn_i_1 <- predict(Churn_Model, Data_i_1, type = 'response')
  
  # Return the output
  return(Churn_i_1)
}

### Churn prediction function for a single period with feedback adjustments
Churn_Pred_Adj_i_1 <- function(P_Old, Churn_Model, Cat_Offer, Data_i_1) {
  # Competitiveness adjustment for previous renewal offers
  Data_i_1$Competitiveness <- Data_i_1$D_1 / P_Old - 1
  Data_i_1$Premium_New_Base <- P_Old * Data_i_1$Risk_Base
  Data_i_1$Undershooting_1 <- pmax(Data_i_1$D_1 - Data_i_1$Premium_New_Base, 
                                   0, na.rm = TRUE)
  Data_i_1$Undershooting_2 <- pmax(Data_i_1$D_2 - Data_i_1$Premium_New_Base, 
                                   0, na.rm = TRUE)
  
  # Offers to consider
  Data_i_1$Offer[] <- Cat_Offer
  
  # Determine churn predictions
  Churn_i_1 <- predict(Churn_Model, Data_i_1, type = 'response')
  
  # Return the output
  return(Churn_i_1)
}

### LP optimization function for a single period
LP_Rsymphony_optim_i_1 <- function(Remain_i_t, Profit_i_t, Churn_i_1, P_Old_i_1, C_i_1, RHS_con, 
                                   Sign_con, Sparse_con, alpha, Offers) {
  # Preliminary variables
  T_max <- length(Offers)
  N_1 <- ncol(Churn_i_1)
  
  # Objective coefficients of LP, where we can eliminate one decision variable per renewal
  # since the sum of the decision variables needs to equal one for each renewal
  LP_Profit <- Remain_i_t * t(1 - Churn_i_1) * 
    (P_Old_i_1 %o% (1 + Offers) - C_i_1) + t(Profit_i_t)
  LP_obj <- as.vector(t(LP_Profit[, -T_max] - LP_Profit[, T_max]))
  
  # Constraint coefficients and signs of LP
  LP_RHS <- Churn_i_1[T_max, ]
  LP_con_sparse <- rbind(Sparse_con[1:N_1, ], as.vector(t(t(Churn_i_1[-T_max, ]) - LP_RHS)) / N_1)
  LP_RHS <- c(head(RHS_con, -1), alpha - sum(LP_RHS) / N_1)
  LP_sign <- Sign_con
  
  # Optimization of LP
  LP_opt <- Rsymphony_solve_LP(max = TRUE, obj = LP_obj, dir = LP_sign, rhs = LP_RHS, 
                               mat = LP_con_sparse, types = rep('B', length(LP_obj)),
                               verbosity = -2)
  
  # Optimal offers
  LP_Ind <- sapply(1:N_1, function(i) min(T_max, which(LP_opt$solution[
    (1 + (T_max - 1) * (i - 1)):((T_max - 1) * i)] == 1)))
  LP_Offers <- Offers[LP_Ind]
  
  # Optimal churn and profit
  LP_Churn <- Churn_i_1[cbind(LP_Ind, 1:N_1)]
  LP_Profit <- LP_Profit[cbind(1:N_1, LP_Ind)]
  
  # Return the output
  return(rbind(LP_opt$status, LP_Ind, LP_Offers, LP_Churn, LP_Profit))
}

### NLP constrained objective function for all periods
NLP_obj_con_i_tau <- function(Offers, P_Old, Offer_Grid, Data_i_tau, Ind_t, Churn_Model, 
                              Costs_i, eps, alpha) {
  # Categorical offers to consider
  Cat_Offers <- as.integer(round(Offers))
  
  # Premium old
  P_Old_2 <- P_Old * (1 + Offer_Grid[Cat_Offers[Ind_t[[1]]]])
  P_Old_3 <- P_Old_2 * (1 + Offer_Grid[Cat_Offers[Ind_t[[2]]]])
  
  # Estimated churn rate
  Churn_i_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                              Cat_Offer = as.character(Cat_Offers[Ind_t[[1]]]), 
                              Data_i_1 = Data_i_tau[Ind_t[[1]], ])
  Churn_i_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_2, Churn_Model = Logit_GLM, 
                                  Cat_Offer = as.character(Cat_Offers[Ind_t[[2]]]), 
                                  Data_i_1 = Data_i_tau[Ind_t[[2]], ])
  Churn_i_3 <- Churn_Pred_Adj_i_1(P_Old = P_Old_3, Churn_Model = Logit_GLM, 
                                  Cat_Offer = as.character(Cat_Offers[Ind_t[[3]]]), 
                                  Data_i_1 = Data_i_tau[Ind_t[[3]], ])
  
  # Expected profit
  Profit_i_1 <- (1 - Churn_i_1) * (P_Old_2 - Costs_i[[1]])
  Profit_i_2 <- (1 - Churn_i_1) * (1 - Churn_i_2) * (P_Old_3 - Costs_i[[2]])
  Profit_i_3 <- (1 - Churn_i_1) * (1 - Churn_i_2) * (1 - Churn_i_3) * 
    (P_Old_3 * (1 + Offer_Grid[Cat_Offers[Ind_t[[3]]]]) - Costs_i[[3]])
  
  # Objective value of NLP
  NLP_obj <- -sum(Profit_i_1 + Profit_i_2 + Profit_i_3)
  
  # Check if the inequality constraints are met
  NLP_con <- c(sum(Churn_i_1), sum(Churn_i_2), sum(Churn_i_3)) / length(P_Old) - alpha
  NLP_obj_con <- (sum(NLP_con <= eps) == length(alpha)) * NLP_obj
  
  # Return the output
  return(NLP_obj_con)
}

### Multi-year horizon
# Customer and period indices
Periods_Max <- 3
Indiv <- unique(AMI$Index_Cust)
N_Indiv <- length(Indiv)
# Filter on customers with at least Periods_Max renewals
Filter_N_i <- sapply(1:N_Indiv, function(i) sum(AMI$Index_Cust == Indiv[i]))
Indiv_Adj <- Indiv[which(Filter_N_i >= Periods_Max)]
N_Indiv_Adj <- length(Indiv_Adj)
AMI_tau <- AMI[which(AMI$Index_Cust %in% Indiv_Adj), ]
AMI_i <- AMI_tau[1, ]
AMI_i[] <- NULL
for (i in 1:N_Indiv_Adj) {
  AMI_Temp <- AMI_tau[which(AMI_tau$Index_Cust == Indiv_Adj[i]), ]
  AMI_i <- rbind(AMI_i, tail(AMI_Temp[
    order(AMI_Temp$Expiry_Date, decreasing = FALSE), ], Periods_Max))
}

### Period inputs
# Indicators
Ind_t_1 <- 1 + Periods_Max * (1:N_Indiv_Adj - 1)
Ind_t_2 <- 2 + Periods_Max * (1:N_Indiv_Adj - 1)
Ind_t_3 <- 3 + Periods_Max * (1:N_Indiv_Adj - 1)
# Base premium
P_Old_Base <- AMI_i$Premium_Old[Ind_t_1]
# Costs
Costs_1 <- AMI_i$Costs[Ind_t_1]
Costs_2 <- AMI_i$Costs[Ind_t_2]
Costs_3 <- AMI_i$Costs[Ind_t_3]
# Data
Data_i_1 <- AMI_i[Ind_t_1, ]
Data_i_2 <- AMI_i[Ind_t_2, ]
Data_i_3 <- AMI_i[Ind_t_3, ]

### Churn predictions - Period 1
Offer_Grid <- as.numeric(Q_means)
Churn_i_1 <- t(sapply(1:T_max, function(t_1)
  Churn_Pred_i_1(Churn_Model = Logit_GLM, Cat_Offer = as.character(t_1), Data_i_1 = Data_i_1)))
Profit_i_1 <- (1 - Churn_i_1) * t(P_Old_Base %o% (1 + Offer_Grid) - Costs_1)
alpha_Range_1 <- apply(sapply(1:N_Indiv_Adj, function(i) 
  c(min(Churn_i_1[, i]), Churn_i_1[which.max(Profit_i_1[, i]), i])), 1, mean)

### Churn of actual offers
Churn_Act_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                              Cat_Offer = as.character(AMI_i$Offer[Ind_t_1]), 
                              Data_i_1 = Data_i_1)
Churn_Act_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[as.numeric(as.character(AMI_i$Offer[Ind_t_1]))]), 
                                  Churn_Model = Logit_GLM, 
                                  Cat_Offer = as.character(AMI_i$Offer[Ind_t_2]),
                                  Data_i_1 = Data_i_2)
Churn_Act_3 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[as.numeric(as.character(AMI_i$Offer[Ind_t_1]))]) *
                                    (1 + Offer_Grid[as.numeric(as.character(AMI_i$Offer[Ind_t_2]))]), 
                                  Churn_Model = Logit_GLM, 
                                  Cat_Offer = as.character(AMI_i$Offer[Ind_t_3]),
                                  Data_i_1 = Data_i_3)
Churn_Act_mean <- c(mean(Churn_Act_1), mean(Churn_Act_2), mean(Churn_Act_3))

### LP setup
# Seed for the random number generator
seed <- 20200101
# Sparse matrix of nrow constraints and ncol variables
LP_con_sparse <- Matrix(0, nrow = N_Indiv_Adj + 1, ncol = N_Indiv_Adj * (T_max - 1), 
                        byrow = TRUE, sparse = TRUE)
LP_con_sparse[cbind(rep(1:N_Indiv_Adj, each = T_max - 1), 1:(N_Indiv_Adj * (T_max - 1)))] <- 1
LP_con_sparse[N_Indiv_Adj + 1, ] <- 1
# Constraints' signs
LP_sign <- c(rep('<=', N_Indiv_Adj), '<=')
# Maximal admissible overall portfolio churn rate
alpha_t <- Churn_Act_mean
# Constraints' RHSs
LP_RHS <- c(rep(1, N_Indiv_Adj), NA)
# Constraint tolerance
LP_eps <- 1e-6

##### LP multi-period optimization #####
### Pre-allocate output - Period 2
Status_2 <- matrix(NA, nrow = T_max, ncol = 4)
Results_Offers_2 <- array(NA, dim = c(T_max, 3, N_Indiv_Adj))

### Period 2
for (t_1 in 1:T_max) {
  # Premium old - Period 2
  P_Old_2 <- P_Old_Base * (1 + Offer_Grid[t_1])
  # Churn - Period 2
  Churn_i_2 <- t(sapply(1:T_max, function(t_2) 
    Churn_Pred_Adj_i_1(P_Old = P_Old_2, Churn_Model = Logit_GLM, 
                       Cat_Offer = as.character(t_2), Data_i_1 = Data_i_2)))
  # Profit - Period 2
  Profit_i_2 <- (1 - Churn_i_1[t_1, ]) * (1 - Churn_i_2) * 
    t(P_Old_2 %o% (1 + Offer_Grid) - Costs_2)
  # Pre-allocate output - Period 3
  Status_3 <- matrix(NA, nrow = T_max, ncol = 2)
  Results_Offers_3 <- array(NA, dim = c(T_max, 2, N_Indiv_Adj))
  
  ### Period 3
  for (t_2 in 1:T_max) {
    # Premium old - Period 3
    P_Old_3 <- P_Old_2 * (1 + Offer_Grid[t_2])
    # Churn - Period 3
    Churn_i_3 <- t(sapply(1:T_max, function(t_3) 
      Churn_Pred_Adj_i_1(P_Old = P_Old_3, Churn_Model = Logit_GLM, 
                         Cat_Offer = as.character(t_3), Data_i_1 = Data_i_3)))
    # Ensure reproducibility by setting the seed for the random number generator
    set.seed(seed)
    # LP optimization - Period 3
    LP_opt_3 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1[t_1, ]) * (1 - Churn_i_2[t_2, ]), 
                                       Profit_i_t = matrix(0, nrow = T_max, ncol = N_Indiv_Adj), 
                                       Churn_i_1 = Churn_i_3, P_Old_i_1 = P_Old_3, 
                                       C_i_1 = Costs_3, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[3], 
                                       Offers = Offer_Grid)
    # Store results if feasible - Period 3
    Status_3[t_2, ] <- c(LP_opt_3[1, 1], mean(LP_opt_3[4, ]))
    if (Status_3[t_2, 2] <= (alpha_t[3] + LP_eps)) {
      Results_Offers_3[t_2, , ] <- rbind(LP_opt_3[3, ], Profit_i_1[t_1, ] + Profit_i_2[t_2, ] +
                                           LP_opt_3[5, ])
    }
  }
  
  # Ensure reproducibility by setting the seed for the random number generator
  set.seed(seed)
  # LP optimization - Period 2
  LP_opt_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1[t_1, ]), 
                                     Profit_i_t = Results_Offers_3[, 2, ], 
                                     Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                     C_i_1 = Costs_2, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                     Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                     Offers = Offer_Grid) 
  # Store results if feasible - Period 2
  Status_2[t_1, ] <- c(LP_opt_2[1, 1], mean(LP_opt_2[4, ]), 
                       apply(Status_3[LP_opt_2[2, ], ], 2, max, na.rm = TRUE))
  if (Status_2[t_1, 2] <= (alpha_t[2] + LP_eps)) {
    Results_Offers_2[t_1, , ] <- rbind(LP_opt_2[3, ], 
                                       Results_Offers_3[cbind(LP_opt_2[2, ], 1, 1:N_Indiv_Adj)],
                                       Results_Offers_3[cbind(LP_opt_2[2, ], 2, 1:N_Indiv_Adj)])
  }
}

# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# LP optimization - Period 1
LP_opt_1 <- LP_Rsymphony_optim_i_1(Remain_i_t = 1, 
                                   Profit_i_t = Results_Offers_2[, 3, ], 
                                   Churn_i_1 = Churn_i_1, P_Old_i_1 = P_Old_Base, 
                                   C_i_1 = Costs_1, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                   Sparse_con = LP_con_sparse, alpha = alpha_t[1], 
                                   Offers = Offer_Grid) 
# Store results if feasible - Period 1
Status_1 <- c(LP_opt_1[1, 1], mean(LP_opt_1[4, ]), 
              apply(Status_2[LP_opt_1[2, ], ], 2, max, na.rm = TRUE))
if (Status_1[2] <= (alpha_t[1] + LP_eps)) {
  Results_Offers_1 <- rbind(LP_opt_1[3, ], 
                            Results_Offers_2[cbind(LP_opt_1[2, ], 1, 1:N_Indiv_Adj)],
                            Results_Offers_2[cbind(LP_opt_1[2, ], 2, 1:N_Indiv_Adj)],
                            Results_Offers_2[cbind(LP_opt_1[2, ], 3, 1:N_Indiv_Adj)])
}

### LP implied churn and profit
# Indices
Ind_LP_1 <- sapply(1:ncol(Results_Offers_1), function(i) which(Offer_Grid == Results_Offers_1[1, i]))
Ind_LP_2 <- sapply(1:ncol(Results_Offers_1), function(i) which(Offer_Grid == Results_Offers_1[2, i]))
Ind_LP_3 <- sapply(1:ncol(Results_Offers_1), function(i) which(Offer_Grid == Results_Offers_1[3, i]))
# Churn
Churn_LP_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, Cat_Offer = as.character(Ind_LP_1), 
                             Data_i_1 = Data_i_1)
Churn_LP_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Results_Offers_1[1, ]), 
                                 Churn_Model = Logit_GLM, Cat_Offer = as.character(Ind_LP_2),
                                 Data_i_1 = Data_i_2)
Churn_LP_3 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Results_Offers_1[1, ]) *
                                   (1 + Results_Offers_1[2, ]), 
                                 Churn_Model = Logit_GLM, Cat_Offer = as.character(Ind_LP_3),
                                 Data_i_1 = Data_i_3)
# Profit
Profit_LP <- list()
Profit_LP[[1]] <- (1 - Churn_LP_1) * (P_Old_Base * (1 + Results_Offers_1[1, ]) - Costs_1)
Profit_LP[[2]] <- (1 - Churn_LP_1) * (1 - Churn_LP_2) * (P_Old_Base * (1 + Results_Offers_1[1, ]) *
                                                           (1 + Results_Offers_1[2, ]) - Costs_2)
Profit_LP[[3]] <- (1 - Churn_LP_1) * (1 - Churn_LP_2) * (1 - Churn_LP_3) * 
  (P_Old_Base * (1 + Results_Offers_1[1, ]) * (1 + Results_Offers_1[2, ]) * (1 + Results_Offers_1[3, ]) - Costs_3)

##### NLP multi-period optimization #####
### Preliminary variables
Ind_t <- list(Ind_t_1, Ind_t_2, Ind_t_3)
Costs_i <- list(Costs_1, Costs_2, Costs_3)
NLP_eps <- 1e-6

### NLP initial values of previous optimization
NLP_Init <- rep(NA, nrow(Results_Offers_1) * Periods_Max)
NLP_Init[Ind_t_1] <- Results_Offers_1[1, ]
NLP_Init[Ind_t_2] <- Results_Offers_1[2, ]
NLP_Init[Ind_t_3] <- Results_Offers_1[3, ]
NLP_Ind <- as.numeric(sapply(1:length(NLP_Init), function(i) which(Offer_Grid == NLP_Init[i])))

### NLP optimization algorithm
NLP_Opts <- list('algorithm' = 'NLOPT_LN_SBPLX', 'ftol_rel' = 0, 'ftol_abs' = 0,
                 'xtol_rel' = 0, 'xtol_abs' = rep(0, length(NLP_Init)),
                 'maxeval' = 1000000, 'print_level' = 0)

### NLP lower and upper bounds of decision variables
NLP_LB <- rep(1, length(NLP_Init))
NLP_UB <- rep(as.numeric(T_max), length(NLP_Init))

### NLP objective function with inequality constraints
NLP_obj_con <- function(x) {
  NLP_obj_con_i_tau(Offers = x, P_Old = P_Old_Base, Offer_Grid = Offer_Grid, 
                    Data_i_tau = AMI_i, Ind_t = Ind_t, Churn_Model = Logit_GLM, 
                    Costs_i = Costs_i, eps = NLP_eps, alpha = alpha_t)
}

### NLP optimization
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Retrieve NLP solution
NLP_opt <- nloptr(x0 = NLP_Ind, eval_f = NLP_obj_con,
                  lb = NLP_LB, ub = NLP_UB, opts = NLP_Opts)

### NLP results
NLP_sol <- as.integer(round(NLP_opt$solution))
NLP_profit <- NLP_opt$objective
NLP_status <- NLP_opt$status

### NLP implied churn and profit
# Churn
Churn_NLP_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                              Cat_Offer = as.character(NLP_sol[Ind_t_1]), 
                              Data_i_1 = Data_i_1)
Churn_NLP_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]), 
                                  Churn_Model = Logit_GLM, 
                                  Cat_Offer = as.character(NLP_sol[Ind_t_2]),
                                  Data_i_1 = Data_i_2)
Churn_NLP_3 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]) *
                                    (1 + Offer_Grid[NLP_sol[Ind_t_2]]), 
                                  Churn_Model = Logit_GLM, 
                                  Cat_Offer = as.character(NLP_sol[Ind_t_3]),
                                  Data_i_1 = Data_i_3)
# Profit
Profit_NLP <- list()
Profit_NLP[[1]] <- (1 - Churn_NLP_1) * (P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]) - Costs_1)
Profit_NLP[[2]] <- (1 - Churn_NLP_1) * (1 - Churn_NLP_2) * (P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]) *
                                                              (1 + Offer_Grid[NLP_sol[Ind_t_2]]) - Costs_2)
Profit_NLP[[3]] <- (1 - Churn_NLP_1) * (1 - Churn_NLP_2) * (1 - Churn_NLP_3) * 
  (P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]) * (1 + Offer_Grid[NLP_sol[Ind_t_2]]) * (1 + Offer_Grid[NLP_sol[Ind_t_3]]) - Costs_3)

### Note: To make sure the best possible solution to the multi-period optimization problem
###       has been found, perform an additional LP multi-period optimization

##### Additional LP multi-period optimization #####
### LP initial values of previous optimization
Churn_i_1 <- Churn_NLP_1
Churn_i_2 <- Churn_NLP_2
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]) *
                       (1 + Offer_Grid[NLP_sol[Ind_t_2]]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_3)))

### Re-iterate over the periods
# Period 3
P_Old_3 <- P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]) * (1 + Offer_Grid[NLP_sol[Ind_t_2]])
LP_extra_3 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), 
                                     Profit_i_t = matrix(0, nrow = T_max, ncol = N_Indiv_Adj), 
                                     Churn_i_1 = Churn_i_3, P_Old_i_1 = P_Old_3, C_i_1 = Costs_3, 
                                     RHS_con = LP_RHS, Sign_con = LP_sign, 
                                     Sparse_con = LP_con_sparse, alpha = alpha_t[3], 
                                     Offers = Offer_Grid)
# Period 2
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]]) * (1 + Offer_Grid[t]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
P_Old_2 <- P_Old_Base * (1 + Offer_Grid[NLP_sol[Ind_t_1]])
LP_extra_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1), 
                                     Profit_i_t = t((1 - Churn_i_1) * t((1 - Churn_i_2) * (1 - Churn_i_3)) * 
                                                      (P_Old_2 %o% (1 + Offer_Grid) * (1 + LP_extra_3[3, ]) - Costs_3)), 
                                     Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                     C_i_1 = Costs_2,  RHS_con = LP_RHS, Sign_con = LP_sign, 
                                     Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                     Offers = Offer_Grid)
# Period 1
Churn_i_1 <- t(sapply(1:T_max, function(t)
  Churn_Pred_i_1(Churn_Model = Logit_GLM, Cat_Offer = as.character(t), Data_i_1 = Data_i_1)))
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]), Churn_Model = Logit_GLM, 
                     Cat_Offer = as.character(LP_extra_2[2, ]), Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]) * (1 + LP_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
LP_extra_1 <- LP_Rsymphony_optim_i_1(Remain_i_t = 1, 
                                     Profit_i_t = t(t(1 - Churn_i_1) * (t(1 - Churn_i_2) * 
                                                                          (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_extra_2[3, ]) - Costs_2) +
                                                                          t(1 - Churn_i_3) * 
                                                                          (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_extra_2[3, ]) * (1 + LP_extra_3[3, ]) - Costs_3))),
                                     Churn_i_1 = Churn_i_1, P_Old_i_1 = P_Old_Base, 
                                     C_i_1 = Costs_1, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                     Sparse_con = LP_con_sparse, alpha = alpha_t[1], 
                                     Offers = Offer_Grid)
# Period 2
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_extra_1[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_extra_1[3, ]) * (1 + Offer_Grid[t]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
P_Old_2 <- P_Old_Base * (1 + LP_extra_1[3, ])
LP_2_extra_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_Loop_1), 
                                       Profit_i_t = t((1 - Churn_Loop_1) * t((1 - Churn_i_2) * (1 - Churn_i_3)) * 
                                                        (P_Old_2 %o% (1 + Offer_Grid) * (1 + LP_extra_3[3, ]) - Costs_3)),
                                       Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                       C_i_1 = Costs_2, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                       Offers = Offer_Grid)
# Period 3
Churn_i_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_extra_1[3, ]), 
                                Churn_Model = Logit_GLM, 
                                Cat_Offer = as.character(LP_2_extra_2[2, ]), 
                                Data_i_1 = Data_i_2)
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_extra_1[3, ]) * (1 + LP_2_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), Data_i_1 = Data_i_3)))
P_Old_3 <- P_Old_Base * (1 + LP_extra_1[3, ]) * (1 + LP_2_extra_2[3, ])
LP_2_extra_3 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_Loop_1) * (1 - Churn_i_2),
                                       Profit_i_t = matrix(0, nrow = T_max, ncol = N_Indiv_Adj), 
                                       Churn_i_1 = Churn_i_3, P_Old_i_1 = P_Old_3, 
                                       C_i_1 = Costs_3, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[3], 
                                       Offers = Offer_Grid) 
# Period 2
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_extra_1[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_extra_1[3, ]) * (1 + Offer_Grid[t]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_2_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
P_Old_2 <- P_Old_Base * (1 + LP_extra_1[3, ])
LP_3_extra_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_Loop_1), 
                                       Profit_i_t = t((1 - Churn_Loop_1) * t((1 - Churn_i_2) * (1 - Churn_i_3)) * 
                                                        (P_Old_2 %o% (1 + Offer_Grid) * (1 + LP_2_extra_3[3, ]) - Costs_3)),
                                       Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                       C_i_1 = Costs_2, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                       Offers = Offer_Grid)
# Period 1
Churn_i_1 <- t(sapply(1:T_max, function(t)
  Churn_Pred_i_1(Churn_Model = Logit_GLM, Cat_Offer = as.character(t), Data_i_1 = Data_i_1)))
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]), Churn_Model = Logit_GLM, 
                     Cat_Offer = as.character(LP_3_extra_2[2, ]), Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]) * (1 + LP_3_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_2_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
LP_2_extra_1 <- LP_Rsymphony_optim_i_1(Remain_i_t = 1, 
                                       Profit_i_t = t(t(1 - Churn_i_1) * (t(1 - Churn_i_2) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_3_extra_2[3, ]) - Costs_2) +
                                                                            t(1 - Churn_i_3) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_3_extra_2[3, ]) * (1 + LP_2_extra_3[3, ]) - Costs_3))),
                                       Churn_i_1 = Churn_i_1, P_Old_i_1 = P_Old_Base, 
                                       C_i_1 = Costs_1, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[1], 
                                       Offers = Offer_Grid)
# Period 2
Churn_i_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                            Cat_Offer = as.character(LP_2_extra_1[2, ]), 
                            Data_i_1 = Data_i_1)
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_2_extra_1[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_2_extra_1[3, ]) * (1 + Offer_Grid[t]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_2_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
P_Old_2 <- P_Old_Base * (1 + LP_2_extra_1[3, ])
LP_4_extra_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1), 
                                       Profit_i_t = t((1 - Churn_i_1) * t((1 - Churn_i_2) * (1 - Churn_i_3)) * 
                                                        (P_Old_2 %o% (1 + Offer_Grid) * (1 + LP_2_extra_3[3, ]) - Costs_3)),
                                       Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                       C_i_1 = Costs_2, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                       Offers = Offer_Grid)
# Period 3
Churn_i_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                            Cat_Offer = as.character(LP_2_extra_1[2, ]), 
                            Data_i_1 = Data_i_1)
Churn_i_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_2_extra_1[3, ]), 
                                Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_4_extra_2[2, ]), 
                                Data_i_1 = Data_i_2)
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_2_extra_1[3, ]) * (1 + LP_4_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_3)))
P_Old_3 <- P_Old_Base * (1 + LP_2_extra_1[3, ]) * (1 + LP_4_extra_2[3, ])
LP_3_extra_3 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), 
                                       Profit_i_t = matrix(0, nrow = T_max, ncol = N_Indiv_Adj), 
                                       Churn_i_1 = Churn_i_3, P_Old_i_1 = P_Old_3, 
                                       C_i_1 = Costs_3, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[3], 
                                       Offers = Offer_Grid) 
# Period 2
Churn_i_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                            Cat_Offer = as.character(LP_2_extra_1[2, ]), 
                            Data_i_1 = Data_i_1)
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_2_extra_1[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_2_extra_1[3, ]) * (1 + Offer_Grid[t]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_3_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
P_Old_2 <- P_Old_Base * (1 + LP_2_extra_1[3, ])
LP_5_extra_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1), 
                                       Profit_i_t = t((1 - Churn_i_1) * t((1 - Churn_i_2) * (1 - Churn_i_3)) * 
                                                        (P_Old_2 %o% (1 + Offer_Grid) * (1 + LP_3_extra_3[3, ]) - Costs_3)),
                                       Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                       C_i_1 = Costs_2, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                       Offers = Offer_Grid)
# Period 1
Churn_i_1 <- t(sapply(1:T_max, function(t)
  Churn_Pred_i_1(Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                 Data_i_1 = Data_i_1)))
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]), Churn_Model = Logit_GLM, 
                     Cat_Offer = as.character(LP_5_extra_2[2, ]), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]) * (1 + LP_5_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_3_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
LP_3_extra_1 <- LP_Rsymphony_optim_i_1(Remain_i_t = 1, 
                                       Profit_i_t = t(t(1 - Churn_i_1) * (t(1 - Churn_i_2) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_5_extra_2[3, ]) - Costs_2) +
                                                                            t(1 - Churn_i_3) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_5_extra_2[3, ]) * (1 + LP_3_extra_3[3, ]) - Costs_3))),
                                       Churn_i_1 = Churn_i_1, P_Old_i_1 = P_Old_Base, 
                                       C_i_1 = Costs_1, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[1], 
                                       Offers = Offer_Grid)
# Period 1
Churn_i_1 <- t(sapply(1:T_max, function(t)
  Churn_Pred_i_1(Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                 Data_i_1 = Data_i_1)))
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]), Churn_Model = Logit_GLM, 
                     Cat_Offer = as.character(LP_5_extra_2[2, ]), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]) * (1 + LP_5_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_3_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
LP_4_extra_1 <- LP_Rsymphony_optim_i_1(Remain_i_t = 1, 
                                       Profit_i_t = t(t(1 - Churn_i_1) * (t(1 - Churn_i_2) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_5_extra_2[3, ]) - Costs_2) +
                                                                            t(1 - Churn_i_3) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_5_extra_2[3, ]) * (1 + LP_3_extra_3[3, ]) - Costs_3))),
                                       Churn_i_1 = Churn_i_1, P_Old_i_1 = P_Old_Base, 
                                       C_i_1 = Costs_1, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[1], 
                                       Offers = Offer_Grid)
# Period 2
Churn_i_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                            Cat_Offer = as.character(LP_3_extra_1[2, ]), 
                            Data_i_1 = Data_i_1)
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_3_extra_1[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_3_extra_1[3, ]) * (1 + Offer_Grid[t]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_3_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
P_Old_2 <- P_Old_Base * (1 + LP_3_extra_1[3, ])
LP_6_extra_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1), 
                                       Profit_i_t = t((1 - Churn_i_1) * t((1 - Churn_i_2) * (1 - Churn_i_3)) * 
                                                        (P_Old_2 %o% (1 + Offer_Grid) * (1 + LP_3_extra_3[3, ]) - Costs_3)),
                                       Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                       C_i_1 = Costs_2, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                       Offers = Offer_Grid)
# Period 3
Churn_i_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                            Cat_Offer = as.character(LP_3_extra_1[2, ]), 
                            Data_i_1 = Data_i_1)
Churn_i_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_3_extra_1[3, ]), 
                                Churn_Model = Logit_GLM, 
                                Cat_Offer = as.character(LP_6_extra_2[2, ]), 
                                Data_i_1 = Data_i_2)
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_3_extra_1[3, ]) * (1 + LP_6_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_3)))
P_Old_3 <- P_Old_Base * (1 + LP_3_extra_1[3, ]) * (1 + LP_6_extra_2[3, ])
LP_4_extra_3 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), 
                                       Profit_i_t = matrix(0, nrow = T_max, ncol = N_Indiv_Adj), 
                                       Churn_i_1 = Churn_i_3, P_Old_i_1 = P_Old_3, 
                                       C_i_1 = Costs_3, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[3], 
                                       Offers = Offer_Grid) 
# Period 2
Churn_i_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                            Cat_Offer = as.character(LP_3_extra_1[2, ]), 
                            Data_i_1 = Data_i_1)
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_3_extra_1[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_3_extra_1[3, ]) * (1 + Offer_Grid[t]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_4_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
P_Old_2 <- P_Old_Base * (1 + LP_3_extra_1[3, ])
LP_7_extra_2 <- LP_Rsymphony_optim_i_1(Remain_i_t = (1 - Churn_i_1), 
                                       Profit_i_t = t((1 - Churn_i_1) * t((1 - Churn_i_2) * (1 - Churn_i_3)) * 
                                                        (P_Old_2 %o% (1 + Offer_Grid) * (1 + LP_4_extra_3[3, ]) - Costs_3)),
                                       Churn_i_1 = Churn_i_2, P_Old_i_1 = P_Old_2, 
                                       C_i_1 = Costs_2, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[2], 
                                       Offers = Offer_Grid)
# Period 1
Churn_i_1 <- t(sapply(1:T_max, function(t)
  Churn_Pred_i_1(Churn_Model = Logit_GLM, Cat_Offer = as.character(t), 
                 Data_i_1 = Data_i_1)))
Churn_i_2 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]), Churn_Model = Logit_GLM, 
                     Cat_Offer = as.character(LP_7_extra_2[2, ]), 
                     Data_i_1 = Data_i_2)))
Churn_i_3 <- t(sapply(1:T_max, function(t)
  Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + Offer_Grid[t]) * (1 + LP_7_extra_2[3, ]), 
                     Churn_Model = Logit_GLM, Cat_Offer = as.character(LP_4_extra_3[2, ]), 
                     Data_i_1 = Data_i_3)))
LP_4_extra_1 <- LP_Rsymphony_optim_i_1(Remain_i_t = 1, 
                                       Profit_i_t = t(t(1 - Churn_i_1) * (t(1 - Churn_i_2) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_7_extra_2[3, ]) - Costs_2) +
                                                                            t(1 - Churn_i_3) * 
                                                                            (P_Old_Base %o% (1 + Offer_Grid) * (1 + LP_7_extra_2[3, ]) * (1 + LP_4_extra_3[3, ]) - Costs_3))),
                                       Churn_i_1 = Churn_i_1, P_Old_i_1 = P_Old_Base, 
                                       C_i_1 = Costs_1, RHS_con = LP_RHS, Sign_con = LP_sign, 
                                       Sparse_con = LP_con_sparse, alpha = alpha_t[1], 
                                       Offers = Offer_Grid)

### Additional LP implied churn and profit
# Churn
Churn_LP_Add_1 <- Churn_Pred_i_1(Churn_Model = Logit_GLM, 
                                 Cat_Offer = as.character(LP_4_extra_1[2, ]), 
                                 Data_i_1 = Data_i_1)
Churn_LP_Add_2 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_4_extra_1[3, ]), 
                                     Churn_Model = Logit_GLM, 
                                     Cat_Offer = as.character(LP_7_extra_2[2, ]),
                                     Data_i_1 = Data_i_2)
Churn_LP_Add_3 <- Churn_Pred_Adj_i_1(P_Old = P_Old_Base * (1 + LP_4_extra_1[3, ]) *
                                       (1 + LP_7_extra_2[3, ]), 
                                     Churn_Model = Logit_GLM, 
                                     Cat_Offer = as.character(LP_4_extra_3[2, ]),
                                     Data_i_1 = Data_i_3)
# Profit
Profit_LP_Add <- list()
Profit_LP_Add[[1]] <- (1 - Churn_LP_Add_1) * (P_Old_Base * (1 + LP_4_extra_1[3, ]) - Costs_1)
Profit_LP_Add[[2]] <- (1 - Churn_LP_Add_1) * (1 - Churn_LP_Add_2) * 
  (P_Old_Base * (1 + LP_4_extra_1[3, ]) * (1 + LP_7_extra_2[3, ]) - Costs_2)
Profit_LP_Add[[3]] <- (1 - Churn_LP_Add_1) * (1 - Churn_LP_Add_2) * (1 - Churn_LP_Add_3) * 
  (P_Old_Base * (1 + LP_4_extra_1[3, ]) * (1 + LP_7_extra_2[3, ]) * (1 + LP_4_extra_3[3, ]) - Costs_3)

##### Results #####
### Comparison of optimization offers
cbind(c('LP', 'NLP', 'Additional LP'), 
      c(sum(Profit_LP[[1]]) + sum(Profit_LP[[2]]) + sum(Profit_LP[[3]]),
        sum(Profit_NLP[[1]]) + sum(Profit_NLP[[2]]) + sum(Profit_NLP[[3]]),
        sum(Profit_LP_Add[[1]]) + sum(Profit_LP_Add[[2]]) + sum(Profit_LP_Add[[3]])))
rbind(c('alpha', 'LP', 'NLP', 'Additional LP'),
      cbind(alpha_t, 
            c(mean(Churn_LP_1), mean(Churn_LP_2), mean(Churn_LP_3)), 
            c(mean(Churn_NLP_1), mean(Churn_NLP_2), mean(Churn_NLP_3)),
            c(mean(Churn_LP_Add_1), mean(Churn_LP_Add_2), mean(Churn_LP_Add_3))))

### Tabulate optimal and actual offers
# Pre-allocate output
Opt_Offer <- matrix(rep(NA, N_Indiv_Adj * 2 * Periods_Max), 
                    nrow = N_Indiv_Adj, ncol = 2 * Periods_Max)
# Actual offers
Opt_Offer[, c(Periods_Max + 1:Periods_Max)] <- matrix(AMI_i$Rate_Change,
                                                      nrow = N_Indiv_Adj,
                                                      ncol = Periods_Max,
                                                      byrow = TRUE)
# Optimal offers
Opt_Offer[, 1:Periods_Max] <- cbind(LP_4_extra_1[3, ], LP_7_extra_2[3, ], LP_4_extra_3[3, ])

### Optimal and actual churn
# Pre-allocate output
Opt_Churn <- matrix(rep(NA, N_Indiv_Adj * 2 * Periods_Max), 
                    nrow = N_Indiv_Adj, ncol = 2 * Periods_Max)
# Actual churn
Opt_Churn[, c(Periods_Max + 1:Periods_Max)] <- cbind(Churn_Act_1, Churn_Act_2, Churn_Act_3)
# Optimal churn
Opt_Churn[, 1:Periods_Max] <- cbind(Churn_LP_Add_1, Churn_LP_Add_2, Churn_LP_Add_3)

### Store the final results
