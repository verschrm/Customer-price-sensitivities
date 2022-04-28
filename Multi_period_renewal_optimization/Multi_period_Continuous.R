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
### Load output from conditional dose-response function and efficient frontier procedure
# Load results from Frontier_Continuous.R

##### (N)LP settings #####
### NLP churn estimates for a single period without feedback adjustments
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

### NLP churn estimates for a single period with feedback adjustments
NLP_Churn_i_t <- function(P_Old_t, GPS_Model, Offers, sd_GPS, Data_i_t, Churn_Model) {
  # Competitiveness adjustment for previous renewal offers
  Data_i_t$Competitiveness <- Data_i_t$D_1 / P_Old_t - 1
  Data_i_t$Premium_New_Base <- P_Old_t * Data_i_t$Risk_Base
  Data_i_t$Undershooting_1 <- pmax(Data_i_t$D_1 - Data_i_t$Premium_New_Base, 
                                   0, na.rm = TRUE)
  Data_i_t$Undershooting_2 <- pmax(Data_i_t$D_2 - Data_i_t$Premium_New_Base, 
                                   0, na.rm = TRUE)
  
  # Generalized propensity score input
  Var_t_XGB <- as.matrix(dummy_cols(Data_i_t[, c('Policy_Type', 'Risk_Group', 'Undershooting_1',
                                                 'Undershooting_2', 'Premium_New_Base')], 
                                    remove_most_frequent_dummy = FALSE)[, -c(1, 2)])
  Var_t_XGB <- Matrix(Var_t_XGB, sparse = TRUE)
  Mean_t <- predict(GPS_XGB, newdata = as.matrix(Var_t_XGB), type = 'response')
  Normalizer_t <- apply(sapply(1:(length(Q_bins) - 1), function(t) 
    pnorm((Q_bins[t + 1] - Mean_t) / sd_GPS[t], mean = 0, sd = 1) - 
      pnorm((Q_bins[t] - Mean_t) / sd_GPS[t], mean = 0, sd = 1)), 1, sum)
  
  # Offers to consider
  Data_i_t$Offer <- Offers
  sd_t <- sd_GPS[as.numeric(cut(Offers, breaks = Q_bins, include.lowest = TRUE))]
  
  # Estimated GPS
  Data_i_t$GPS <- dnorm((Offers - Mean_t) / sd_t, mean = 0, sd = 1) / Normalizer_t / sd_t
  
  # Determine churn predictions
  Var_t_XGB <- as.matrix(dummy_cols(
    Data_i_t[, c('Competitiveness', 'Offer', 'GPS', 'Policy_Type', 'Risk_Group',
                 'Undershooting_1', 'Undershooting_2', 'Premium_New_Base')],
    remove_most_frequent_dummy = FALSE)[, -c(4, 5)])
  Var_t_XGB <- Matrix(Var_t_XGB, sparse = TRUE)
  Churn_i_t <- predict(Churn_Model, newdata = as.matrix(Var_t_XGB), type = 'response')
  
  # Return the output
  return(Churn_i_t)
}

### LP local optimization function for a single period
LP_L_Rsymphony_optim_i_t <- function(Offer_i_t, P_Old_i_1, Data_i_t, GPS_Model,
                                     Remain_i_t, Data_i_1, sd_GPS, Mean, Normalizer, 
                                     Churn_Model, RHS_con, Sign_con, Sparse_con, alpha, 
                                     Offers_Init, Num, Scale, LB, UB, t) {
  # Preliminary variables
  N_1 <- length(Offers_Init)
  Scale_Num <- Scale / (Num - 1)
  
  # Form local search domain
  Scale_L_temp <- ifelse((Offers_Init - Scale / 2) < LB, Offers_Init - LB, Scale / 2)
  Scale_U_temp <- ifelse((Offers_Init + Scale / 2) > UB, UB - Offers_Init, Scale / 2)
  Scale_L <- Scale_L_temp + (Scale / 2 - Scale_U_temp)
  Scale_U <- Scale_U_temp + (Scale / 2 - Scale_L_temp)
  Offers <- mapply(function(t, l, u) c(t - l, seq(from = t - ceiling(l / Scale_Num - 1) * Scale_Num, 
                                                  to = t + floor(u / Scale_Num - 1) * Scale_Num,
                                                  length.out = Num - 2), t + u), 
                   Offers_Init, Scale_L, Scale_U)
  
  # Determine corresponding churn predictions
  if (t == 1) {
    Churn_i_1 <- t(sapply(1:Num, function(t) 
      NLP_Churn_i_1(Offers = Offers[t, ], sd_GPS = sd_GPS, Data_i_1 = Data_i_1, 
                    Mean = Mean, Normalizer = Normalizer, Churn_Model = Churn_Model)))
  } else {
    Churn_i_1 <- t(sapply(1:Num, function(t) 
      NLP_Churn_i_t(P_Old_t = P_Old_i_1, GPS_Model = GPS_Model, Offers = Offers[t, ], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_1, Churn_Model = Churn_Model)))
  }
  
  # Adjust profit for feedback effects
  if (t == 3) {
    Profit_i_t <- matrix(0, nrow = N_1, ncol = Num)
  } else if (t == 2) {
    P_t_3 <- P_Old_i_1 * (1 + t(Offers))
    Churn_t_3 <- t(sapply(1:Num, function(t) 
      NLP_Churn_i_t(P_Old_t = P_t_3[, t], GPS_Model = GPS_Model, Offers = Offer_i_t[[1]], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_t[[1]], Churn_Model = Churn_Model)))
    Profit_i_t <- Remain_i_t * t((1 - Churn_i_1) * (1 - Churn_t_3)) * 
      (P_t_3 * (1 + Offer_i_t[[1]]) - Data_i_t[[1]]$Costs)
  } else {
    P_t_2 <- P_Old_i_1 * (1 + t(Offers))
    P_t_3 <- P_t_2 * (1 + Offer_i_t[[2]])
    Churn_t_2 <- t(sapply(1:Num, function(t) 
      NLP_Churn_i_t(P_Old_t = P_t_2[, t], GPS_Model = GPS_Model, Offers = Offer_i_t[[2]], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_t[[2]], Churn_Model = Churn_Model)))
    Churn_t_3 <- t(sapply(1:Num, function(t) 
      NLP_Churn_i_t(P_Old_t = P_t_3[, t], GPS_Model = GPS_Model, Offers = Offer_i_t[[1]], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_t[[1]], Churn_Model = Churn_Model)))
    Profit_i_t <- Remain_i_t * t((1 - Churn_i_1) * (1 - Churn_t_2)) * ((P_t_3 - Data_i_t[[2]]$Costs) +
        t(1 - Churn_t_3) * (P_t_3 * (1 + Offer_i_t[[1]]) - Data_i_t[[1]]$Costs))
  }
  
  # Objective coefficients of LP, where we can eliminate one decision variable per renewal
  # since the sum of the decision variables needs to equal one for each renewal
  LP_Profit <- Remain_i_t * t(1 - Churn_i_1) * (P_Old_i_1 * t(1 + Offers) - 
                                                  Data_i_1$Costs) + Profit_i_t
  LP_obj <- -as.vector(t(LP_Profit[, -Num] - LP_Profit[, Num]))
  
  # Constraint coefficients and signs of LP
  LP_RHS <- Churn_i_1[Num, ]
  LP_con_sparse <- rbind(Sparse_con[1:N_1, ], as.vector(t(t(Churn_i_1[-Num, ]) - LP_RHS)) / N_1)
  LP_RHS <- c(head(RHS_con, -1), alpha - sum(LP_RHS) / N_1)
  LP_sign <- Sign_con
  
  # Optimization of LP
  LP_opt <- Rsymphony_solve_LP(max = FALSE, obj = LP_obj, dir = LP_sign, rhs = LP_RHS, 
                               mat = LP_con_sparse, types = rep('B', length(LP_obj)),
                               time_limit = 30, verbosity = -2)
  
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

### GPS adjustments
GPS_Adj_i_t <- function(P_Old_t, GPS_Model, sd_GPS, Data_i_t) {
  # Competitiveness adjustment for previous renewal offers
  Data_i_t$Competitiveness <- Data_i_t$D_1 / P_Old_t - 1
  Data_i_t$Premium_New_Base <- P_Old_t * Data_i_t$Risk_Base
  Data_i_t$Undershooting_1 <- pmax(Data_i_t$D_1 - Data_i_t$Premium_New_Base, 
                                   0, na.rm = TRUE)
  Data_i_t$Undershooting_2 <- pmax(Data_i_t$D_2 - Data_i_t$Premium_New_Base, 
                                   0, na.rm = TRUE)
  
  # Generalized propensity score input
  Var_t_XGB <- as.matrix(dummy_cols(Data_i_t[, c('Policy_Type', 'Risk_Group', 'Undershooting_1',
                                                 'Undershooting_2', 'Premium_New_Base')], 
                                    remove_most_frequent_dummy = FALSE)[, -c(1, 2)])
  Var_t_XGB <- Matrix(Var_t_XGB, sparse = TRUE)
  Mean_t <- predict(GPS_XGB, newdata = as.matrix(Var_t_XGB), type = 'response')
  Normalizer_t <- apply(sapply(1:(length(Q_bins) - 1), function(t) 
    pnorm((Q_bins[t + 1] - Mean_t) / sd_GPS[t], mean = 0, sd = 1) - 
      pnorm((Q_bins[t] - Mean_t) / sd_GPS[t], mean = 0, sd = 1)), 1, sum)
  
  # Return the output
  return(cbind(Mean_t, Normalizer_t))
}

### LP optimization wrapper
LP_opt_wrapper <- function(Offer_i_t, P_Old_i_1, Data_i_t, GPS_Model,
                           Remain_i_t, Data_i_1, sd_GPS, Mean, Normalizer, 
                           Churn_Model, RHS_con, Sign_con, Sparse_con, alpha, 
                           Offers_Init, Profit_Init, Num, Scales, LB, UB, t, seed, eps) {
  # Preliminary variables
  LP_Num <- length(Scales)
  Init_n <- Offers_Init
  Profit_n <- Profit_Init
  Status_n <- -1
  
  # Loop over local optimizations
  for (n in 1:LP_Num) {
    # Ensure reproducibility by setting the seed for the random number generator
    set.seed(seed)
    # Run the n-th LP optimization
    LP_opt_n <- LP_L_Rsymphony_optim_i_t(Offer_i_t = Offer_i_t, P_Old_i_1 = P_Old_i_1, 
                                         Data_i_t = Data_i_t, GPS_Model = GPS_Model,
                                         Remain_i_t = Remain_i_t, Data_i_1 = Data_i_1, 
                                         sd_GPS = sd_GPS, Mean = Mean, 
                                         Normalizer = Normalizer, Churn_Model = Churn_Model, 
                                         RHS_con = RHS_con, Sign_con = LP_sign, 
                                         Sparse_con = Sparse_con, alpha = alpha, 
                                         Offers_Init = Init_n, Num = Num, Scale = Scales[n], 
                                         LB = LB, UB = UB, t = t)
    
    # If solution has improved and is feasible, update initial values
    if (((-sum(LP_opt_n[4, ])) < Profit_n) & 
        ((floor(mean(LP_opt_n[3, ]) * 1e6) / 1e6) <= (alpha + eps))) {
      Status_n <- LP_opt_n[1, 1]
      Profit_n <- -sum(LP_opt_n[4, ])
      Init_n <- LP_opt_n[2, ]
    }
  }
  
  # Determine the corresponding churn predictions
  if (t == 1) {
    Churn_i_n <- NLP_Churn_i_1(Offers = Init_n, sd_GPS = sd_GPS, Data_i_1 = Data_i_1, 
                               Mean = Mean, Normalizer = Normalizer, Churn_Model = Churn_Model)
  } else {
    Churn_i_n <- NLP_Churn_i_t(P_Old_t = P_Old_i_1, GPS_Model = GPS_Model, Offers = Init_n, 
                               sd_GPS = sd_GPS, Data_i_t = Data_i_1, Churn_Model = Churn_Model)
  }
  
  # Adjust profit for feedback effects
  if (t == 3) {
    Profit_i_t <- rep(0, length(Init_n))
  } else if (t == 2) {
    P_n_3 <- P_Old_i_1 * (1 + Init_n)
    Churn_n_3 <- NLP_Churn_i_t(P_Old_t = P_n_3, GPS_Model = GPS_Model, Offers = Offer_i_t[[1]], 
                               sd_GPS = sd_GPS, Data_i_t = Data_i_t[[1]], Churn_Model = Churn_Model)
    Profit_i_t <- Remain_i_t * (1 - Churn_i_n) * (1 - Churn_n_3) * (P_n_3 * (1 + Offer_i_t[[1]]) - 
                                                                      Data_i_t[[1]]$Costs)
  } else {
    P_n_2 <- P_Old_i_1 * (1 + Init_n)
    P_n_3 <- P_n_2 * (1 + Offer_i_t[[2]])
    Churn_n_2 <- NLP_Churn_i_t(P_Old_t = P_n_2, GPS_Model = GPS_Model, Offers = Offer_i_t[[2]], 
                               sd_GPS = sd_GPS, Data_i_t = Data_i_t[[2]], Churn_Model = Churn_Model)
    Churn_n_3 <- NLP_Churn_i_t(P_Old_t = P_n_3, GPS_Model = GPS_Model, Offers = Offer_i_t[[1]], 
                               sd_GPS = sd_GPS, Data_i_t = Data_i_t[[1]], Churn_Model = Churn_Model)
    Profit_i_t <- (1 - Churn_i_n) * (1 - Churn_n_2) * ((P_n_3 - Data_i_t[[2]]$Costs) +
                                                         (1 - Churn_n_3) * (P_n_3 * (1 + Offer_i_t[[1]]) - Data_i_t[[1]]$Costs))
  }
  Profit_i_n <- Remain_i_t * (1 - Churn_i_n) * (P_Old_i_1 * (1 + Init_n) - 
                                                   Data_i_1$Costs) + Profit_i_t
  
  # Return the output
  return(rbind(Status_n, Init_n, Churn_i_n, Profit_i_n))
}

### Initial local search over offers
ILS_grid_search <- function(ILS_grid, Step_Shift, ILS_Ind, Step_Num, sd_GPS, Data_i_1,
                            Mean, Normalizer, Churn_Model, P_1, C_1, 
                            Offer_i_t, Data_i_t, GPS_Model, Remain_i_t, t) {
  # Preliminary variables
  ILS_Min <- min(ILS_grid) + max(ILS_Ind, 0) * Step_Shift
  ILS_Max <- max(ILS_grid) + min(ILS_Ind, 0) * Step_Shift
  Step_Num_Adj <- ceiling(Step_Num / (max(ILS_grid) - min(ILS_grid)) * (ILS_Max - ILS_Min))
  
  # Initial local search offers
  Offers_ILS <- seq(from = ILS_Min, to = ILS_Max, length.out = Step_Num_Adj)
  
  # Determine the corresponding churn predictions
  if (t == 1) {
    Churn_ILS <- t(sapply(1:Step_Num_Adj, function(t) 
      NLP_Churn_i_1(Offers = Offers_ILS[t], sd_GPS = sd_GPS, Data_i_1 = Data_i_1, 
                    Mean = Mean, Normalizer = Normalizer, Churn_Model = Churn_Model)))
  } else {
    Churn_ILS <- t(sapply(1:Step_Num_Adj, function(t) 
      NLP_Churn_i_t(P_Old_t = P_1, GPS_Model = GPS_Model, Offers = Offers_ILS[t], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_1, Churn_Model = Churn_Model)))
  }
  
  # Adjust profit for feedback effects
  if (t == 3) {
    Profit_i_t <- matrix(0, nrow = length(P_1), ncol = Step_Num_Adj)
  } else if (t == 2) {
    P_3 <- P_1 %o% (1 + Offers_ILS)
    Churn_3 <- t(sapply(1:Step_Num_Adj, function(t) 
      NLP_Churn_i_t(P_Old_t = P_3[, t], GPS_Model = GPS_Model, Offers = Offer_i_t[[1]], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_t[[1]], Churn_Model = Churn_Model)))
    Profit_i_t <- Remain_i_t * t((1 - Churn_ILS) * (1 - Churn_3)) * (P_3 * (1 + Offer_i_t[[1]]) - 
                                                                       Data_i_t[[1]]$Costs)
  } else {
    P_2 <- P_1 %o% (1 + Offers_ILS)
    P_3 <- P_2 * (1 + Offer_i_t[[2]])
    Churn_2 <- t(sapply(1:Step_Num_Adj, function(t) 
      NLP_Churn_i_t(P_Old_t = P_2[, t], GPS_Model = GPS_Model, Offers = Offer_i_t[[2]], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_t[[2]], Churn_Model = Churn_Model)))
    Churn_3 <- t(sapply(1:Step_Num_Adj, function(t) 
      NLP_Churn_i_t(P_Old_t = P_3[, t], GPS_Model = GPS_Model, Offers = Offer_i_t[[1]], 
                    sd_GPS = sd_GPS, Data_i_t = Data_i_t[[1]], Churn_Model = Churn_Model)))
    Profit_i_t <- t((1 - Churn_ILS) * (1 - Churn_2)) * ((P_3 - Data_i_t[[2]]$Costs) +
      t(1 - Churn_3) * (P_3 * (1 + Offer_i_t[[1]]) - Data_i_t[[1]]$Costs))
  }
  Profit_ILS <- t(Remain_i_t * t(1 - Churn_ILS) * (P_1 %o% (1 + Offers_ILS) - C_1) + Profit_i_t)
  
  # Minimal and maximal indices
  Ind_ILS <- sapply(1:nrow(Data_i_1), function(i) 
    c(which.min(Churn_ILS[, i]), which.max(Profit_ILS[, i])))
  
  # Return the output
  return(list(cbind(Offers_ILS[Ind_ILS[1, ]], Offers_ILS[Ind_ILS[2, ]]),
              cbind(rbind(mean(Churn_ILS[cbind(Ind_ILS[1, ], 1:nrow(Data_i_1))]), 
                          -sum(Profit_ILS[cbind(Ind_ILS[1, ], 1:nrow(Data_i_1))])),
              rbind(mean(Churn_ILS[cbind(Ind_ILS[2, ], 1:nrow(Data_i_1))]), 
                    -sum(Profit_ILS[cbind(Ind_ILS[2, ], 1:nrow(Data_i_1))])))))
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

### Generalized propensity score input - Period 1
Var_1_XGB <- as.matrix(dummy_cols(Data_i_1[, c('Policy_Type', 'Risk_Group', 'Undershooting_1',
                                               'Undershooting_2', 'Premium_New_Base')], 
                                  remove_most_frequent_dummy = FALSE)[, -c(1, 2)])
Var_1_XGB <- Matrix(Var_1_XGB, sparse = TRUE)
Mean_1 <- predict(GPS_XGB, newdata = as.matrix(Var_1_XGB), type = 'response')
Normalizer_1 <- apply(sapply(1:(length(Q_bins) - 1), function(t) 
  pnorm((Q_bins[t + 1] - Mean_1) / sd_GPS_XGB[t], mean = 0, sd = 1) - 
    pnorm((Q_bins[t] - Mean_1) / sd_GPS_XGB[t], mean = 0, sd = 1)), 1, sum)

### Churn of actual offers
Churn_Act_1 <- NLP_Churn_i_1(Offers = AMI_i$Offer[Ind_t_1], sd_GPS = sd_GPS_XGB, 
                             Data_i_1 = Data_i_1, Mean = Mean_1, Normalizer = Normalizer_1, 
                             Churn_Model = Churn_XGB)
Churn_Act_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + AMI_i$Offer[Ind_t_1]), 
                             GPS_Model = GPS_XGB, Offers = AMI_i$Offer[Ind_t_2],
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2,
                             Churn_Model = Churn_XGB)
Churn_Act_3 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + AMI_i$Offer[Ind_t_1])  * 
                               (1 + AMI_i$Offer[Ind_t_2]),
                             GPS_Model = GPS_XGB, Offers = AMI_i$Offer[Ind_t_3],
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3,
                             Churn_Model = Churn_XGB)
Churn_Act_mean <- c(mean(Churn_Act_1), mean(Churn_Act_2), mean(Churn_Act_3))

### (N)LP setup
# Seed for the random number generator
seed <- 20200101
# LP local search domain properties
Num_L <- 3
Scales_L <- 2 * c(5, 4, 3, 2, 1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 
                  0.02, 0.01, 0.005, 0.004, 0.003, 0.002, 0.001) / 100
# LP sparse matrix of nrow constraints and ncol variables
LP_con_sparse_L <- Matrix(0, nrow = N_Indiv_Adj, ncol = N_Indiv_Adj * (Num_L - 1), 
                          byrow = TRUE, sparse = TRUE)
LP_con_sparse_L[cbind(rep(1:N_Indiv_Adj, each = Num_L - 1), 1:(N_Indiv_Adj * (Num_L - 1)))] <- 1
LP_con_sparse_L <- rbind(LP_con_sparse_L, 1)
# Constraints' signs
LP_sign <- c(rep('<=', N_Indiv_Adj), '<=')
# Maximal admissible overall portfolio churn rate
alpha_t <- Churn_Act_mean
# Constraints' RHSs
LP_RHS <- c(rep(1, N_Indiv_Adj), NA)
# (N)LP constraint tolerance
LP_eps <- 0
NLP_eps <- 1e-6

# NLP global and local lower and upper bounds of decision variables
NLP_LB_G <- rep(as.numeric(Q_bins[1]), nrow(AMI_i))
NLP_UB_G <- rep(as.numeric(tail(Q_bins, 1)), nrow(AMI_i))
NLP_LB_L <- rep(as.numeric(Q_bins[1]), N_Indiv_Adj)
NLP_UB_L <- rep(as.numeric(tail(Q_bins, 1)), N_Indiv_Adj)

### Initial local search options
Step_Shift <- 2 / 100
Tries <- ceiling((max(T_Offers) - min(T_Offers) - Step_Shift) / Step_Shift)
Step_Num <- ceiling(5 * 100 * (max(T_Offers) - min(T_Offers)))

### Initial churn grid
Step_Grid <- 1 / 100
Offer_Grid <- seq(from = ceiling(min(T_Offers) / Step_Grid) * Step_Grid, 
                  to = floor(max(T_Offers) / Step_Grid) * Step_Grid, by = Step_Grid)
Offer_Grid <- c(min(T_Offers), Offer_Grid, max(T_Offers))
T_t <- length(Offer_Grid)
sd_Grid <- sd_GPS_XGB[as.numeric(cut(Offer_Grid, breaks = Q_bins, include.lowest = TRUE))]

##### LP multi-period optimization #####
### Pre-allocate output - Period 2
Status_2 <- matrix(NA, nrow = T_t, ncol = 6)
Results_Offers_2 <- array(NA, dim = c(T_t, 3, N_Indiv_Adj))

### Churn - Period 1
Churn_i_1 <- t(sapply(1:T_t, function(t_1) 
  NLP_Churn_i_1(Offers = Offer_Grid[t_1], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)))

### Profit - Period 1
Profit_i_1 <- (1 - Churn_i_1) * t(P_Old_Base %o% (1 + Offer_Grid) - Costs_1)

### Period 2
for (t_1 in 1:T_t) {
  # Premium old - Period 2
  P_Old_2 <- P_Old_Base
  P_Old_2 <- P_Old_2 * (1 + Offer_Grid[t_1])
  # Churn - Period 2
  Churn_i_2 <- t(sapply(1:T_t, function(t_2) 
    NLP_Churn_i_t(P_Old_t = P_Old_2, GPS_Model = GPS_XGB, Offers = Offer_Grid[t_2], 
                  sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, Churn_Model = Churn_XGB)))
  # Profit - Period 2
  Profit_i_2 <- (1 - Churn_i_1[t_1, ]) * (1 - Churn_i_2) * 
    t(P_Old_2 %o% (1 + Offer_Grid) - Costs_2)
  # Pre-allocate output - Period 3
  Status_3 <- matrix(NA, nrow = T_t, ncol = 3)
  Results_Offers_3 <- array(NA, dim = c(T_t, 2, N_Indiv_Adj))
  
  ### Period 3
  for (t_2 in 1:T_t) {
    # Premium old - Period 3
    P_Old_3 <- P_Old_2 * (1 + Offer_Grid[t_2])
    # Initial local search - Period 3
    ILS_gs_3 <- sapply((-Tries):Tries, function(ils)
      ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                     Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_3,
                     Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_3, 
                     C_1 = Costs_3, Offer_i_t = NULL, Data_i_t = NULL, GPS_Model = GPS_XGB, 
                     Remain_i_t = (1 - Churn_i_1[t_1, ]) * (1 - Churn_i_2[t_2, ]), t = 3))
    # Summarize initial local search results - Period 3
    Summary_gs_3 <- matrix(unlist(ILS_gs_3[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
    ILS_gs_3 <- matrix(unlist(ILS_gs_3[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
    # Select closest, optimal (and feasible) initial solution - Period 3
    Iter_gs_3 <- which(Summary_gs_3[1, ] <= (alpha_t[3] + NLP_eps))
    Ind_gs_3 <- Iter_gs_3[which.min(Summary_gs_3[2, Iter_gs_3])]
    # Initial values - Period 3
    LP_Init_3 <- ILS_gs_3[, Ind_gs_3]
    LP_Profit_3 <- Summary_gs_3[2, Ind_gs_3]
    # LP optimizations - Period 3
    LP_opt_3 <- LP_opt_wrapper(Offer_i_t = NULL, P_Old_i_1 = P_Old_3, Data_i_t = NULL, 
                               GPS_Model = GPS_XGB, 
                               Remain_i_t = (1 - Churn_i_1[t_1, ]) * (1 - Churn_i_2[t_2, ]), 
                               Data_i_1 = Data_i_3, sd_GPS = sd_GPS_XGB, Mean = NULL, 
                               Normalizer = NULL, Churn_Model = Churn_XGB, 
                               RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[3], 
                               Offers_Init = LP_Init_3, Profit_Init = LP_Profit_3, 
                               Num = Num_L, Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 3, seed = seed, 
                               eps = LP_eps)
    # Store results if feasible and an improvement - Period 3
    if ((mean(LP_opt_3[3, ]) <= (alpha_t[3] + NLP_eps)) & ((-sum(LP_opt_3[4, ])) < LP_Profit_3)) {
      Status_3[t_2, ] <- c(LP_opt_3[1, 1], mean(LP_opt_3[3, ]), -sum(LP_opt_3[4, ]))
      Results_Offers_3[t_2, , ] <- rbind(LP_opt_3[2, ], 
                                         Profit_i_1[t_1, ] + Profit_i_2[t_2, ] + LP_opt_3[4, ])
    # Otherwise stick with the initial values - Period 3
    } else {
      Churn_i_3 <- NLP_Churn_i_t(P_Old_t = P_Old_3, GPS_Model = GPS_XGB, Offers = LP_Init_3, 
                                 sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, Churn_Model = Churn_XGB)
      LP_Profit_3 <- (1 - Churn_i_1[t_1, ]) * (1 - Churn_i_2[t_2, ]) * (1 - Churn_i_3) * 
        (P_Old_3 * (1 + LP_Init_3) - Costs_3)
      Status_3[t_2, ] <- c(-1, mean(Churn_i_3), -sum(LP_Profit_3))
      Results_Offers_3[t_2, , ] <- rbind(LP_Init_3, 
                                         Profit_i_1[t_1, ] + Profit_i_2[t_2, ] + LP_Profit_3)
    }
  }
  
  # Best initial starting point - Period 2
  LP_Ind_2 <- which.min(-apply(Results_Offers_3[, 2, ], 1, sum))
  # Initial local grid search - Period 2
  ILS_gs_2 <- sapply((-Tries):Tries, function(ils)
    ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                   Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                   Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                   C_1 = Costs_2, Offer_i_t = list(Results_Offers_3[LP_Ind_2, 1, ]), 
                   Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                   Remain_i_t = (1 - Churn_i_1[t_1, ]), t = 2))
  # Summarize initial local search results - Period 2
  Summary_gs_2 <- matrix(unlist(ILS_gs_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
  ILS_gs_2 <- matrix(unlist(ILS_gs_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
  # Select closest, optimal (and feasible) initial solution - Period 2
  Iter_gs_2 <- which(Summary_gs_2[1, ] <= (alpha_t[2] + NLP_eps))
  Ind_gs_2 <- Iter_gs_2[which.min(Summary_gs_2[2, Iter_gs_2])]
  # Initial values - Period 2
  LP_Init_2 <- ILS_gs_2[, Ind_gs_2]
  LP_Profit_2 <- Summary_gs_2[2, Ind_gs_2]
  # LP optimizations - Period 2
  LP_opt_2 <- LP_opt_wrapper(Offer_i_t = list(Results_Offers_3[LP_Ind_2, 1, ]), 
                             P_Old_i_1 = P_Old_2, Data_i_t = list(Data_i_3), 
                             GPS_Model = GPS_XGB, Remain_i_t = (1 - Churn_i_1[t_1, ]), 
                             Data_i_1 = Data_i_2, sd_GPS = sd_GPS_XGB, Mean = NULL, 
                             Normalizer = NULL, Churn_Model = Churn_XGB, 
                             RHS_con = LP_RHS, Sign_con = LP_sign, 
                             Sparse_con = LP_con_sparse_L, alpha = alpha_t[2], 
                             Offers_Init = LP_Init_2, Profit_Init = LP_Profit_2, 
                             Num = Num_L, Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                             UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                             eps = LP_eps)
  # Store results if feasible and an improvement - Period 2
  if ((mean(LP_opt_2[3, ]) <= (alpha_t[2] + NLP_eps)) & 
      ((-sum(LP_opt_2[4, ])) < LP_Profit_2)) {
    Status_2[t_1, ] <- c(LP_opt_2[1, 1], mean(LP_opt_2[3, ]), -sum(LP_opt_2[4, ]), 
                         Status_3[LP_Ind_2, ])
    Results_Offers_2[t_1, , ] <- rbind(LP_opt_2[2, ], Results_Offers_3[LP_Ind_2, 1, ],
                                       Profit_i_1[t_1, ] + LP_opt_2[4, ])
  # Otherwise stick with the initial values - Period 2
  } else {
    Churn_i_2 <- NLP_Churn_i_t(P_Old_t = P_Old_2, GPS_Model = GPS_XGB, Offers = LP_Init_2, 
                               sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, 
                               Churn_Model = Churn_XGB)
    Churn_i_3 <- NLP_Churn_i_t(P_Old_t = P_Old_2 * (1 + LP_Init_2), GPS_Model = GPS_XGB, 
                               Offers = Results_Offers_3[LP_Ind_2, 1, ], 
                               sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, 
                               Churn_Model = Churn_XGB)
    LP_Profit_2 <- (1 - Churn_i_1[t_1, ]) * (1 - Churn_i_2) * 
      ((P_Old_2 * (1 + LP_Init_2) - Costs_2) + (1 - Churn_i_3) * 
         (P_Old_2 * (1 + LP_Init_2) * (1 + Results_Offers_3[LP_Ind_2, 1, ]) - Costs_3))
    Status_2[t_1, ] <- c(-1, mean(Churn_i_2), -sum(LP_Profit_2), Status_3[LP_Ind_2, ])
    Results_Offers_2[t_1, , ] <- rbind(LP_Init_2, Results_Offers_3[LP_Ind_2, 1, ],
                                       Profit_i_1[t_1, ] + LP_Profit_2)
  }
}

# Best initial starting point - Period 1
LP_Ind_1 <- which.min(-apply(Results_Offers_2[, 3, ], 1, sum))
# Initial local search - Period 1
ILS_gs_1 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1,
                 Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB, 
                 P_1 = P_Old_Base, C_1 = Costs_1, 
                 Offer_i_t = list(Results_Offers_2[LP_Ind_1, 2, ], 
                                  Results_Offers_2[LP_Ind_1, 1, ]), 
                 Data_i_t = list(Data_i_3, Data_i_2), GPS_Model = GPS_XGB, 
                 Remain_i_t = 1, t = 1))
# Summarize initial local search results - Period 1
Summary_gs_1 <- matrix(unlist(ILS_gs_1[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_gs_1 <- matrix(unlist(ILS_gs_1[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Select closest, optimal (and feasible) initial solution - Period 1
Iter_gs_1 <- which(Summary_gs_1[1, ] <= (alpha_t[1] + NLP_eps))
Ind_gs_1 <- Iter_gs_1[which.min(Summary_gs_1[2, Iter_gs_1])]
# Initial values - Period 1
LP_Init_1 <- ILS_gs_1[, Ind_gs_1]
LP_Profit_1 <- Summary_gs_1[2, Ind_gs_1]
# LP optimization - Period 1
LP_opt_1 <- LP_opt_wrapper(Offer_i_t = list(Results_Offers_2[LP_Ind_1, 2, ], 
                                            Results_Offers_2[LP_Ind_1, 1, ]), 
                           P_Old_i_1 = P_Old_Base, Data_i_t = list(Data_i_3, Data_i_2), 
                           GPS_Model = GPS_XGB, Remain_i_t = 1, Data_i_1 = Data_i_1, 
                           sd_GPS = sd_GPS_XGB, Mean = Mean_1, Normalizer = Normalizer_1, 
                           Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                           Sparse_con = LP_con_sparse_L, alpha = alpha_t[1], 
                           Offers_Init = LP_Init_1, Profit_Init = LP_Profit_1, 
                           Num = Num_L, Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                           UB = as.numeric(tail(Q_bins, 1)), t = 1, seed = seed, 
                           eps = LP_eps)
# Store results if feasible and an improvement - Period 1
if ((mean(LP_opt_1[3, ]) <= (alpha_t[1] + NLP_eps)) & 
    ((-sum(LP_opt_1[4, ])) < LP_Profit_1)) {
  Status_1 <- c(LP_opt_1[1, 1], mean(LP_opt_1[3, ]), -sum(LP_opt_1[4, ]), 
                       Status_2[LP_Ind_1, ])
  Results_Offers_1 <- rbind(LP_opt_1[2, ], Results_Offers_2[LP_Ind_1, c(1, 2), ],
                                     LP_opt_1[4, ])
# Otherwise stick with the initial values - Period 2
} else {
  Churn_i_1 <- NLP_Churn_i_1(Offers = LP_Init_1, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                             Mean = Mean_1, Normalizer = Normalizer_1, 
                             Churn_Model = Churn_XGB)
  Churn_i_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + LP_Init_1), GPS_Model = GPS_XGB, 
                             Offers = Results_Offers_2[LP_Ind_1, 1, ], 
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, 
                             Churn_Model = Churn_XGB)
  Churn_i_3 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + LP_Init_1) * 
                               (1 + Results_Offers_2[LP_Ind_1, 1, ]), 
                             GPS_Model = GPS_XGB, Offers = Results_Offers_2[LP_Ind_1, 2, ], 
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, 
                             Churn_Model = Churn_XGB)
  LP_Profit_1 <- (1 - Churn_i_1) * ((P_Old_Base * (1 + LP_Init_1) - Costs_1) +
    (1 - Churn_i_2) * ((P_Old_Base * (1 + LP_Init_1) * 
                          (1 + Results_Offers_2[LP_Ind_1, 1, ]) - Costs_2) +
    (1 - Churn_i_3) * (P_Old_Base * (1 + LP_Init_1) * 
                         (1 + Results_Offers_2[LP_Ind_1, 1, ]) * 
                         (1 + Results_Offers_2[LP_Ind_1, 2, ]) - Costs_3)))
  Status_1 <- c(-1, mean(Churn_i_1), -sum(LP_Profit_1), Status_2[LP_Ind_1, ])
  Results_Offers_1 <- rbind(LP_Init_1, Results_Offers_2[LP_Ind_1, c(1, 2), ],
                                     LP_Profit_1)
}

### LP implied churn and profit
# Churn
Churn_LP_1 <- NLP_Churn_i_1(Offers = Results_Offers_1[1, ], sd_GPS = sd_GPS_XGB, 
                            Data_i_1 = Data_i_1, Mean = Mean_1, Normalizer = Normalizer_1, 
                            Churn_Model = Churn_XGB)
Churn_LP_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + Results_Offers_1[1, ]), 
                            GPS_Model = GPS_XGB, Offers = Results_Offers_1[2, ], 
                            sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, 
                            Churn_Model = Churn_XGB)
Churn_LP_3 <- NLP_Churn_i_t(P_Old = P_Old_Base * (1 + Results_Offers_1[1, ]) * 
                              (1 + Results_Offers_1[2, ]), 
                            GPS_Model = GPS_XGB, Offers = Results_Offers_1[3, ], 
                            sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, 
                            Churn_Model = Churn_XGB)
# Profit
Profit_LP <- list()
Profit_LP[[1]] <- (1 - Churn_LP_1) * (P_Old_Base * (1 + Results_Offers_1[1, ]) - Costs_1)
Profit_LP[[2]] <- (1 - Churn_LP_1) * (1 - Churn_LP_2) * 
  (P_Old_Base * (1 + Results_Offers_1[1, ]) * (1 + Results_Offers_1[2, ]) - Costs_2)
Profit_LP[[3]] <- (1 - Churn_LP_1) * (1 - Churn_LP_2) * (1 - Churn_LP_3) * 
  (P_Old_Base * (1 + Results_Offers_1[1, ]) * (1 + Results_Offers_1[2, ]) * 
     (1 + Results_Offers_1[3, ]) - Costs_3)

##### NLP multi-period optimization #####
### NLP optimization algorithm
NLP_Opts <- list('algorithm' = 'NLOPT_LN_SBPLX', 'ftol_rel' = 0, 'ftol_abs' = 0,
                 'xtol_rel' = 0, 'xtol_abs' = rep(0, nrow(AMI_i)),
                 'maxeval' = 100000, 'print_level' = 0)

### NLP objective function with inequality constraints
NLP_obj_con <- function(x) {
  # Premium old
  P_Old_2 <- P_Old_Base * (1 + x[Ind_t_1])
  P_Old_3 <- P_Old_2 * (1 + x[Ind_t_2])
  
  # Estimated churn rate
  Churn_i_1 <- NLP_Churn_i_1(Offers = x[Ind_t_1], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                             Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
  Churn_i_2 <- NLP_Churn_i_t(P_Old_t = P_Old_2, GPS_Model = GPS_XGB, Offers = x[Ind_t_2], 
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, Churn_Model = Churn_XGB)
  Churn_i_3 <- NLP_Churn_i_t(P_Old_t = P_Old_3, GPS_Model = GPS_XGB, Offers = x[Ind_t_3], 
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, Churn_Model = Churn_XGB)
  
  # Expected profit
  Profit_i_1 <- (1 - Churn_i_1) * (P_Old_2 - Costs_1)
  Profit_i_2 <- (1 - Churn_i_1) * (1 - Churn_i_2) * (P_Old_3 - Costs_2)
  Profit_i_3 <- (1 - Churn_i_1) * (1 - Churn_i_2) * (1 - Churn_i_3) * 
    (P_Old_3 * (1 + x[Ind_t_3]) - Costs_3)
  
  # Objective value of NLP
  NLP_obj <- -sum(Profit_i_1 + Profit_i_2 + Profit_i_3)
  
  # Check if the inequality constraints are met
  NLP_con <- c(sum(Churn_i_1), sum(Churn_i_2), sum(Churn_i_3)) / N_Indiv_Adj
  NLP_obj_con <- (sum(NLP_con <= (alpha_t + NLP_eps)) == length(alpha_t)) * NLP_obj
  
  # Return the output
  return(NLP_obj_con)
}

### NLP initial values of previous optimization
NLP_Init <- rep(NA, nrow(Results_Offers_1) * Periods_Max)
NLP_Init[Ind_t_1] <- Results_Offers_1[1, ]
NLP_Init[Ind_t_2] <- Results_Offers_1[2, ]
NLP_Init[Ind_t_3] <- Results_Offers_1[3, ]

### NLP optimization
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Retrieve NLP solution
NLP_opt <- nloptr(x0 = NLP_Init, eval_f = NLP_obj_con,
                  lb = NLP_LB_G, ub = NLP_UB_G, opts = NLP_Opts)

### NLP results
NLP_sol <- NLP_opt$solution
NLP_profit <- NLP_opt$objective
NLP_status <- NLP_opt$status

### NLP implied churn and profit
# Churn
Churn_NLP_1 <- NLP_Churn_i_1(Offers = NLP_sol[Ind_t_1], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
Churn_NLP_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + NLP_sol[Ind_t_1]), 
                             GPS_Model = GPS_XGB, Offers = NLP_sol[Ind_t_2], 
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, Churn_Model = Churn_XGB)
Churn_NLP_3 <- NLP_Churn_i_t(P_Old = P_Old_Base * (1 + NLP_sol[Ind_t_1]) * (1 + NLP_sol[Ind_t_2]), 
                             GPS_Model = GPS_XGB, Offers = NLP_sol[Ind_t_3], 
                             sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, Churn_Model = Churn_XGB)
# Profit
Profit_NLP <- list()
Profit_NLP[[1]] <- (1 - Churn_NLP_1) * (P_Old_Base * (1 + NLP_sol[Ind_t_1]) - Costs_1)
Profit_NLP[[2]] <- (1 - Churn_NLP_1) * (1 - Churn_NLP_2) * (P_Old_Base * (1 + NLP_sol[Ind_t_1]) *
                                                              (1 + NLP_sol[Ind_t_2]) - Costs_2)
Profit_NLP[[3]] <- (1 - Churn_NLP_1) * (1 - Churn_NLP_2) * (1 - Churn_NLP_3) * 
  (P_Old_Base * (1 + NLP_sol[Ind_t_1]) * (1 + NLP_sol[Ind_t_2]) * (1 + NLP_sol[Ind_t_3]) - Costs_3)

### Note: To make sure the best possible solution to the multi-period optimization problem
###       has been found, perform an additional LP and NLP multi-period optimization

##### Additional LP multi-period optimization #####
### LP initial values of previous optimization
Churn_i_1 <- Churn_NLP_1
Churn_i_2 <- Churn_NLP_2

### Re-iterate over the periods
## Period 3
# Premium old - Period 3
P_Old_3 <- P_Old_Base * (1 + NLP_sol[Ind_t_1]) * (1 + NLP_sol[Ind_t_2])
# Initial local search - Period 3
ILS_gs_extra_3 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_3,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_3, 
                 C_1 = Costs_3, Offer_i_t = NULL, Data_i_t = NULL, GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), t = 3))
Summary_gs_extra_3 <- matrix(unlist(ILS_gs_extra_3[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_gs_extra_3 <- matrix(unlist(ILS_gs_extra_3[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 3
Iter_gs_extra_3 <- which(Summary_gs_extra_3[1, ] <= (alpha_t[3] + NLP_eps))
Ind_gs_extra_3 <- Iter_gs_extra_3[which.min(Summary_gs_extra_3[2, Iter_gs_extra_3])]
LP_Init_extra_3 <- ILS_gs_extra_3[, Ind_gs_extra_3]
LP_Profit_extra_3 <- Summary_gs_extra_3[2, Ind_gs_extra_3]
# LP optimizations - Period 3
LP_extra_3 <- LP_opt_wrapper(Offer_i_t = NULL, P_Old_i_1 = P_Old_3, Data_i_t = NULL, 
                             GPS_Model = GPS_XGB, Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), 
                             Data_i_1 = Data_i_3, sd_GPS = sd_GPS_XGB, Mean = NULL, 
                             Normalizer = NULL, Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                             Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                             alpha = alpha_t[3], Offers_Init = LP_Init_extra_3, 
                             Profit_Init = LP_Profit_extra_3, Num = Num_L, Scales = Scales_L,
                             LB = as.numeric(Q_bins[1]), UB = as.numeric(tail(Q_bins, 1)),
                             t = 3, seed = seed, eps = LP_eps)

## Period 2
# Premium old - Period 2
P_Old_2 <- P_Old_Base * (1 + NLP_sol[Ind_t_1])
# Initial local search - Period 2
ILS_gs_extra_2 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                 C_1 = Costs_2, Offer_i_t = list(LP_extra_3[2, ]), 
                 Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1), t = 2))
Summary_gs_extra_2 <- matrix(unlist(ILS_gs_extra_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_gs_extra_2 <- matrix(unlist(ILS_gs_extra_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 2
Iter_gs_extra_2 <- which(Summary_gs_extra_2[1, ] <= (alpha_t[2] + NLP_eps))
Ind_gs_extra_2 <- Iter_gs_extra_2[which.min(Summary_gs_extra_2[2, Iter_gs_extra_2])]
LP_Init_extra_2 <- ILS_gs_extra_2[, Ind_gs_extra_2]
LP_Profit_extra_2 <- Summary_gs_extra_2[2, Ind_gs_extra_2]
# LP optimizations - Period 2
LP_extra_2 <- LP_opt_wrapper(Offer_i_t = list(LP_extra_3[2, ]), P_Old_i_1 = P_Old_2, 
                             Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB,
                             Remain_i_t = (1 - Churn_i_1), Data_i_1 = Data_i_2, 
                             sd_GPS = sd_GPS_XGB, Mean = NULL, Normalizer = NULL, 
                             Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                             Sparse_con = LP_con_sparse_L, alpha = alpha_t[2], 
                             Offers_Init = LP_Init_extra_2, Profit_Init = LP_Profit_extra_2, 
                             Num = Num_L, Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                             UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                             eps = LP_eps)

## Period 1
# Initial local search - Period 1
ILS_gs_extra_1 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1,
                 Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB, 
                 P_1 = P_Old_Base, C_1 = Costs_1, 
                 Offer_i_t = list(LP_extra_3[2, ], LP_extra_2[2, ]), 
                 Data_i_t = list(Data_i_3, Data_i_2), GPS_Model = GPS_XGB, 
                 Remain_i_t = 1, t = 1))
Summary_gs_extra_1 <- matrix(unlist(ILS_gs_extra_1[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_gs_extra_1 <- matrix(unlist(ILS_gs_extra_1[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 1
Iter_gs_extra_1 <- which(Summary_gs_extra_1[1, ] <= (alpha_t[1] + NLP_eps))
Ind_gs_extra_1 <- Iter_gs_extra_1[which.min(Summary_gs_extra_1[2, Iter_gs_extra_1])]
LP_Init_extra_1 <- ILS_gs_extra_1[, Ind_gs_extra_1]
LP_Profit_extra_1 <- Summary_gs_extra_1[2, Ind_gs_extra_1]
# LP optimizations - Period 1
LP_extra_1 <- LP_opt_wrapper(Offer_i_t = list(LP_extra_3[2, ], LP_extra_2[2, ]), 
                             P_Old_i_1 = P_Old_Base, Data_i_t = list(Data_i_3, Data_i_2), 
                             GPS_Model = GPS_XGB, Remain_i_t = 1, Data_i_1 = Data_i_1, 
                             sd_GPS = sd_GPS_XGB, Mean = Mean_1, Normalizer = Normalizer_1, 
                             Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                             Sparse_con = LP_con_sparse_L, alpha = alpha_t[1], 
                             Offers_Init = LP_Init_extra_1, Profit_Init = LP_Profit_extra_1, 
                             Num = Num_L, Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                             UB = as.numeric(tail(Q_bins, 1)), t = 1, seed = seed, 
                             eps = LP_eps)

## Period 2
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Premium old - Period 2
P_Old_2 <- P_Old_Base * (1 + LP_extra_1[2, ])
# Initial local search - Period 2
ILS_2_gs_extra_2 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                 C_1 = Costs_2, Offer_i_t = list(LP_extra_3[2, ]), 
                 Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1), t = 2))
Summary_2_gs_extra_2 <- matrix(unlist(ILS_2_gs_extra_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_2_gs_extra_2 <- matrix(unlist(ILS_2_gs_extra_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 2
Iter_2_gs_extra_2 <- which(Summary_2_gs_extra_2[1, ] <= (alpha_t[2] + NLP_eps))
Ind_2_gs_extra_2 <- Iter_2_gs_extra_2[which.min(Summary_2_gs_extra_2[2, Iter_2_gs_extra_2])]
LP_2_Init_extra_2 <- ILS_2_gs_extra_2[, Ind_2_gs_extra_2]
LP_2_Profit_extra_2 <- Summary_2_gs_extra_2[2, Ind_2_gs_extra_2]
# LP optimizations - Period 2
LP_2_extra_2 <- LP_opt_wrapper(Offer_i_t = list(LP_extra_3[2, ]), P_Old_i_1 = P_Old_2, 
                               Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB,
                               Remain_i_t = (1 - Churn_i_1), Data_i_1 = Data_i_2, 
                               sd_GPS = sd_GPS_XGB, Mean = NULL, Normalizer = NULL, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[2], 
                               Offers_Init = LP_2_Init_extra_2, 
                               Profit_Init = LP_2_Profit_extra_2, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                               eps = LP_eps)

## Period 3
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Churn - Period 2
Churn_i_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + LP_extra_1[2, ]), 
                           GPS_Model = GPS_XGB, Offers = LP_2_extra_2[2, ], 
                           sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, Churn_Model = Churn_XGB)
# Premium old - Period 3
P_Old_3 <- P_Old_Base * (1 + LP_extra_1[2, ]) * (1 + LP_2_extra_2[2, ])
# Initial local search - Period 3
ILS_2_gs_extra_3 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_3,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_3, 
                 C_1 = Costs_3, Offer_i_t = NULL, Data_i_t = NULL, GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), t = 3))
Summary_2_gs_extra_3 <- matrix(unlist(ILS_2_gs_extra_3[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_2_gs_extra_3 <- matrix(unlist(ILS_2_gs_extra_3[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 3
Iter_2_gs_extra_3 <- which(Summary_2_gs_extra_3[1, ] <= (alpha_t[3] + NLP_eps))
Ind_2_gs_extra_3 <- Iter_2_gs_extra_3[which.min(Summary_2_gs_extra_3[2, Iter_2_gs_extra_3])]
LP_2_Init_extra_3 <- ILS_2_gs_extra_3[, Ind_2_gs_extra_3]
LP_2_Profit_extra_3 <- Summary_2_gs_extra_3[2, Ind_2_gs_extra_3]
# LP optimizations - Period 3
LP_2_extra_3 <- LP_opt_wrapper(Offer_i_t = NULL, P_Old_i_1 = P_Old_3, Data_i_t = NULL, 
                               GPS_Model = GPS_XGB, 
                               Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), 
                               Data_i_1 = Data_i_3, sd_GPS = sd_GPS_XGB, Mean = NULL, 
                               Normalizer = NULL, Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                               Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                               alpha = alpha_t[3], Offers_Init = LP_2_Init_extra_3, 
                               Profit_Init = LP_2_Profit_extra_3, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 3, seed = seed, 
                               eps = LP_eps)

## Period 2
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Premium old - Period 2
P_Old_2 <- P_Old_Base * (1 + LP_extra_1[2, ])
# Initial local search - Period 2
ILS_3_gs_extra_2 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                 C_1 = Costs_2, Offer_i_t = list(LP_2_extra_3[2, ]), 
                 Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1), t = 2))
Summary_3_gs_extra_2 <- matrix(unlist(ILS_3_gs_extra_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_3_gs_extra_2 <- matrix(unlist(ILS_3_gs_extra_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 2
Iter_3_gs_extra_2 <- which(Summary_3_gs_extra_2[1, ] <= (alpha_t[2] + NLP_eps))
Ind_3_gs_extra_2 <- Iter_3_gs_extra_2[which.min(Summary_3_gs_extra_2[2, Iter_3_gs_extra_2])]
LP_3_Init_extra_2 <- ILS_3_gs_extra_2[, Ind_3_gs_extra_2]
LP_3_Profit_extra_2 <- Summary_3_gs_extra_2[2, Ind_3_gs_extra_2]
# LP optimizations - Period 2
LP_3_extra_2 <- LP_opt_wrapper(Offer_i_t = list(LP_2_extra_3[2, ]), P_Old_i_1 = P_Old_2, 
                               Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB,
                               Remain_i_t = (1 - Churn_i_1), Data_i_1 = Data_i_2, 
                               sd_GPS = sd_GPS_XGB, Mean = NULL, Normalizer = NULL, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[2], 
                               Offers_Init = LP_3_Init_extra_2, 
                               Profit_Init = LP_3_Profit_extra_2, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                               eps = LP_eps)

## Period 1
# Initial local search - Period 1
ILS_2_gs_extra_1 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1,
                 Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB, P_1 = P_Old_Base, 
                 C_1 = Costs_1, Offer_i_t = list(LP_2_extra_3[2, ], LP_3_extra_2[2, ]), 
                 Data_i_t = list(Data_i_3, Data_i_2), GPS_Model = GPS_XGB, 
                 Remain_i_t = 1, t = 1))
Summary_2_gs_extra_1 <- matrix(unlist(ILS_2_gs_extra_1[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_2_gs_extra_1 <- matrix(unlist(ILS_2_gs_extra_1[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 1
Iter_2_gs_extra_1 <- which(Summary_2_gs_extra_1[1, ] <= (alpha_t[1] + NLP_eps))
Ind_2_gs_extra_1 <- Iter_2_gs_extra_1[which.min(Summary_2_gs_extra_1[2, Iter_2_gs_extra_1])]
LP_2_Init_extra_1 <- ILS_2_gs_extra_1[, Ind_2_gs_extra_1]
LP_2_Profit_extra_1 <- Summary_2_gs_extra_1[2, Ind_2_gs_extra_1]
# LP optimizations - Period 1
LP_2_extra_1 <- LP_opt_wrapper(Offer_i_t = list(LP_2_extra_3[2, ], LP_3_extra_2[2, ]), 
                               P_Old_i_1 = P_Old_Base, Data_i_t = list(Data_i_3, Data_i_2), 
                               GPS_Model = GPS_XGB, Remain_i_t = 1, Data_i_1 = Data_i_1, 
                               sd_GPS = sd_GPS_XGB, Mean = Mean_1, Normalizer = Normalizer_1, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[1], 
                               Offers_Init = LP_2_Init_extra_1, Profit_Init = LP_2_Profit_extra_1,
                               Num = Num_L, Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 1, seed = seed, 
                               eps = LP_eps)

## Period 2
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_2_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Premium old - Period 2
P_Old_2 <- P_Old_Base * (1 + LP_2_extra_1[2, ])
# Initial local search - Period 2
ILS_4_gs_extra_2 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                 C_1 = Costs_2, Offer_i_t = list(LP_2_extra_3[2, ]), 
                 Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1), t = 2))
Summary_4_gs_extra_2 <- matrix(unlist(ILS_4_gs_extra_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_4_gs_extra_2 <- matrix(unlist(ILS_4_gs_extra_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 2
Iter_4_gs_extra_2 <- which(Summary_4_gs_extra_2[1, ] <= (alpha_t[2] + NLP_eps))
Ind_4_gs_extra_2 <- Iter_4_gs_extra_2[which.min(Summary_4_gs_extra_2[2, Iter_4_gs_extra_2])]
LP_4_Init_extra_2 <- ILS_4_gs_extra_2[, Ind_4_gs_extra_2]
LP_4_Profit_extra_2 <- Summary_4_gs_extra_2[2, Ind_4_gs_extra_2]
# LP optimizations - Period 2
LP_4_extra_2 <- LP_opt_wrapper(Offer_i_t = list(LP_2_extra_3[2, ]), P_Old_i_1 = P_Old_2, 
                               Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB,
                               Remain_i_t = (1 - Churn_i_1), Data_i_1 = Data_i_2, 
                               sd_GPS = sd_GPS_XGB, Mean = NULL, Normalizer = NULL, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[2], 
                               Offers_Init = LP_4_Init_extra_2, Profit_Init = LP_4_Profit_extra_2, 
                               Num = Num_L, Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                               eps = LP_eps)

## Period 3
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_2_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Churn - Period 2
Churn_i_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + LP_2_extra_1[2, ]), 
                           GPS_Model = GPS_XGB, Offers = LP_4_extra_2[2, ], 
                           sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, Churn_Model = Churn_XGB)
# Premium old - Period 3
P_Old_3 <- P_Old_Base * (1 + LP_2_extra_1[2, ]) * (1 + LP_4_extra_2[2, ])
# Initial local search - Period 3
ILS_3_gs_extra_3 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_3,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_3, 
                 C_1 = Costs_3, Offer_i_t = NULL, Data_i_t = NULL, GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), t = 3))
Summary_3_gs_extra_3 <- matrix(unlist(ILS_3_gs_extra_3[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_3_gs_extra_3 <- matrix(unlist(ILS_3_gs_extra_3[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 3
Iter_3_gs_extra_3 <- which(Summary_3_gs_extra_3[1, ] <= (alpha_t[3] + NLP_eps))
Ind_3_gs_extra_3 <- Iter_3_gs_extra_3[which.min(Summary_3_gs_extra_3[2, Iter_3_gs_extra_3])]
LP_3_Init_extra_3 <- ILS_3_gs_extra_3[, Ind_3_gs_extra_3]
LP_3_Profit_extra_3 <- Summary_3_gs_extra_3[2, Ind_3_gs_extra_3]
# LP optimizations - Period 3
LP_3_extra_3 <- LP_opt_wrapper(Offer_i_t = NULL, P_Old_i_1 = P_Old_3, Data_i_t = NULL, 
                               GPS_Model = GPS_XGB, 
                               Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), 
                               Data_i_1 = Data_i_3, sd_GPS = sd_GPS_XGB, Mean = NULL, 
                               Normalizer = NULL, Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                               Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                               alpha = alpha_t[3], Offers_Init = LP_3_Init_extra_3, 
                               Profit_Init = LP_3_Profit_extra_3, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 3, seed = seed, e
                               ps = LP_eps)

## Period 2
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_2_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Premium old - Period 2
P_Old_2 <- P_Old_Base * (1 + LP_2_extra_1[2, ])
# Initial local search - Period 2
ILS_5_gs_extra_2 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                 C_1 = Costs_2, Offer_i_t = list(LP_3_extra_3[2, ]), 
                 Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1), t = 2))
Summary_5_gs_extra_2 <- matrix(unlist(ILS_5_gs_extra_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_5_gs_extra_2 <- matrix(unlist(ILS_5_gs_extra_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 2
Iter_5_gs_extra_2 <- which(Summary_5_gs_extra_2[1, ] <= (alpha_t[2] + NLP_eps))
Ind_5_gs_extra_2 <- Iter_5_gs_extra_2[which.min(Summary_5_gs_extra_2[2, Iter_5_gs_extra_2])]
LP_5_Init_extra_2 <- ILS_5_gs_extra_2[, Ind_5_gs_extra_2]
LP_5_Profit_extra_2 <- Summary_5_gs_extra_2[2, Ind_5_gs_extra_2]
# LP optimizations - Period 2
LP_5_extra_2 <- LP_opt_wrapper(Offer_i_t = list(LP_3_extra_3[2, ]), P_Old_i_1 = P_Old_2, 
                               Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB,
                               Remain_i_t = (1 - Churn_i_1), Data_i_1 = Data_i_2, 
                               sd_GPS = sd_GPS_XGB, Mean = NULL, Normalizer = NULL, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                               Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                               alpha = alpha_t[2], Offers_Init = LP_5_Init_extra_2, 
                               Profit_Init = LP_5_Profit_extra_2, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                               eps = LP_eps)

### Period 1
# Initial local search - Period 1
ILS_3_gs_extra_1 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1,
                 Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB, P_1 = P_Old_Base, 
                 C_1 = Costs_1, Offer_i_t = list(LP_3_extra_3[2, ], LP_5_extra_2[2, ]), 
                 Data_i_t = list(Data_i_3, Data_i_2), GPS_Model = GPS_XGB, 
                 Remain_i_t = 1, t = 1))
Summary_3_gs_extra_1 <- matrix(unlist(ILS_3_gs_extra_1[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_3_gs_extra_1 <- matrix(unlist(ILS_3_gs_extra_1[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 1
Iter_3_gs_extra_1 <- which(Summary_3_gs_extra_1[1, ] <= (alpha_t[1] + NLP_eps))
Ind_3_gs_extra_1 <- Iter_3_gs_extra_1[which.min(Summary_3_gs_extra_1[2, Iter_3_gs_extra_1])]
LP_3_Init_extra_1 <- ILS_3_gs_extra_1[, Ind_3_gs_extra_1]
LP_3_Profit_extra_1 <- Summary_3_gs_extra_1[2, Ind_3_gs_extra_1]
# LP optimizations - Period 1
LP_3_extra_1 <- LP_opt_wrapper(Offer_i_t = list(LP_3_extra_3[2, ], LP_5_extra_2[2, ]), 
                               P_Old_i_1 = P_Old_Base, Data_i_t = list(Data_i_3, Data_i_2), 
                               GPS_Model = GPS_XGB, Remain_i_t = 1, Data_i_1 = Data_i_1, 
                               sd_GPS = sd_GPS_XGB, Mean = Mean_1, Normalizer = Normalizer_1, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[1], 
                               Offers_Init = LP_3_Init_extra_1, 
                               Profit_Init = LP_3_Profit_extra_1, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 1, seed = seed, 
                               eps = LP_eps)

## Period 1
# Initial local search - Period 1
ILS_4_gs_extra_1 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1,
                 Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB, P_1 = P_Old_Base, 
                 C_1 = Costs_1, Offer_i_t = list(LP_3_extra_3[2, ], LP_5_extra_2[2, ]), 
                 Data_i_t = list(Data_i_3, Data_i_2), GPS_Model = GPS_XGB, 
                 Remain_i_t = 1, t = 1))
Summary_4_gs_extra_1 <- matrix(unlist(ILS_4_gs_extra_1[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_4_gs_extra_1 <- matrix(unlist(ILS_4_gs_extra_1[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 1
Iter_4_gs_extra_1 <- which(Summary_4_gs_extra_1[1, ] <= (alpha_t[1] + NLP_eps))
Ind_4_gs_extra_1 <- Iter_4_gs_extra_1[which.min(Summary_4_gs_extra_1[2, Iter_4_gs_extra_1])]
LP_4_Init_extra_1 <- ILS_4_gs_extra_1[, Ind_4_gs_extra_1]
LP_4_Profit_extra_1 <- Summary_4_gs_extra_1[2, Ind_4_gs_extra_1]
# LP optimizations - Period 1
LP_4_extra_1 <- LP_opt_wrapper(Offer_i_t = list(LP_3_extra_3[2, ], LP_5_extra_2[2, ]), 
                               P_Old_i_1 = P_Old_Base, Data_i_t = list(Data_i_3, Data_i_2), 
                               GPS_Model = GPS_XGB, Remain_i_t = 1, Data_i_1 = Data_i_1, 
                               sd_GPS = sd_GPS_XGB, Mean = Mean_1, Normalizer = Normalizer_1, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[1], 
                               Offers_Init = LP_4_Init_extra_1, 
                               Profit_Init = LP_4_Profit_extra_1, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 1, seed = seed, 
                               eps = LP_eps)

## Period 2
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_3_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Premium old - Period 2
P_Old_2 <- P_Old_Base * (1 + LP_3_extra_1[2, ])
# Initial local search - Period 2
ILS_6_gs_extra_2 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                 C_1 = Costs_2, Offer_i_t = list(LP_3_extra_3[2, ]), 
                 Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1), t = 2))
Summary_6_gs_extra_2 <- matrix(unlist(ILS_6_gs_extra_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_6_gs_extra_2 <- matrix(unlist(ILS_6_gs_extra_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 2
Iter_6_gs_extra_2 <- which(Summary_6_gs_extra_2[1, ] <= (alpha_t[2] + NLP_eps))
Ind_6_gs_extra_2 <- Iter_6_gs_extra_2[which.min(Summary_6_gs_extra_2[2, Iter_6_gs_extra_2])]
LP_6_Init_extra_2 <- ILS_6_gs_extra_2[, Ind_6_gs_extra_2]
LP_6_Profit_extra_2 <- Summary_6_gs_extra_2[2, Ind_6_gs_extra_2]
# LP optimizations - Period 2
LP_6_extra_2 <- LP_opt_wrapper(Offer_i_t = list(LP_3_extra_3[2, ]), P_Old_i_1 = P_Old_2, 
                               Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB,
                               Remain_i_t = (1 - Churn_i_1), Data_i_1 = Data_i_2, 
                               sd_GPS = sd_GPS_XGB, Mean = NULL, Normalizer = NULL, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[2], 
                               Offers_Init = LP_6_Init_extra_2, 
                               Profit_Init = LP_6_Profit_extra_2, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                               eps = LP_eps)

## Period 3
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_3_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Churn - Period 2
Churn_i_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + LP_3_extra_1[2, ]), 
                           GPS_Model = GPS_XGB, Offers = LP_6_extra_2[2, ], 
                           sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, Churn_Model = Churn_XGB)
# Premium old - Period 2
P_Old_3 <- P_Old_Base * (1 + LP_3_extra_1[2, ]) * (1 + LP_6_extra_2[2, ])
# Initial local search - Period 3
ILS_4_gs_extra_3 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_3,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_3, 
                 C_1 = Costs_3, Offer_i_t = NULL, Data_i_t = NULL, GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), t = 3))
Summary_4_gs_extra_3 <- matrix(unlist(ILS_4_gs_extra_3[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_4_gs_extra_3 <- matrix(unlist(ILS_4_gs_extra_3[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 3
Iter_4_gs_extra_3 <- which(Summary_4_gs_extra_3[1, ] <= (alpha_t[3] + NLP_eps))
Ind_4_gs_extra_3 <- Iter_4_gs_extra_3[which.min(Summary_4_gs_extra_3[2, Iter_4_gs_extra_3])]
LP_4_Init_extra_3 <- ILS_4_gs_extra_3[, Ind_4_gs_extra_3]
LP_4_Profit_extra_3 <- Summary_4_gs_extra_3[2, Ind_4_gs_extra_3]
# LP optimizations - Period 3
LP_4_extra_3 <- LP_opt_wrapper(Offer_i_t = NULL, P_Old_i_1 = P_Old_3, Data_i_t = NULL, 
                               GPS_Model = GPS_XGB, 
                               Remain_i_t = (1 - Churn_i_1) * (1 - Churn_i_2), 
                               Data_i_1 = Data_i_3, sd_GPS = sd_GPS_XGB, Mean = NULL, 
                               Normalizer = NULL, Churn_Model = Churn_XGB, RHS_con = LP_RHS, 
                               Sign_con = LP_sign, Sparse_con = LP_con_sparse_L, 
                               alpha = alpha_t[3], Offers_Init = LP_4_Init_extra_3, 
                               Profit_Init = LP_4_Profit_extra_3, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 3, seed = seed, 
                               eps = LP_eps)

## Period 2
# Churn - Period 1
Churn_i_1 <- NLP_Churn_i_1(Offers = LP_3_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                           Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
# Premium old - Period 2
P_Old_2 <- P_Old_Base * (1 + LP_3_extra_1[2, ])
# Initial local search - Period 2
ILS_7_gs_extra_2 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_2,
                 Mean = NULL, Normalizer = NULL, Churn_Model = Churn_XGB, P_1 = P_Old_2, 
                 C_1 = Costs_2, Offer_i_t = list(LP_4_extra_3[2, ]), 
                 Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB, 
                 Remain_i_t = (1 - Churn_i_1), t = 2))
Summary_7_gs_extra_2 <- matrix(unlist(ILS_7_gs_extra_2[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_7_gs_extra_2 <- matrix(unlist(ILS_7_gs_extra_2[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 2
Iter_7_gs_extra_2 <- which(Summary_7_gs_extra_2[1, ] <= (alpha_t[2] + NLP_eps))
Ind_7_gs_extra_2 <- Iter_7_gs_extra_2[which.min(Summary_7_gs_extra_2[2, Iter_7_gs_extra_2])]
LP_7_Init_extra_2 <- ILS_7_gs_extra_2[, Ind_7_gs_extra_2]
LP_7_Profit_extra_2 <- Summary_7_gs_extra_2[2, Ind_7_gs_extra_2]
# LP optimizations - Period 2
LP_7_extra_2 <- LP_opt_wrapper(Offer_i_t = list(LP_4_extra_3[2, ]), P_Old_i_1 = P_Old_2, 
                               Data_i_t = list(Data_i_3), GPS_Model = GPS_XGB,
                               Remain_i_t = (1 - Churn_i_1), Data_i_1 = Data_i_2, 
                               sd_GPS = sd_GPS_XGB, Mean = NULL, Normalizer = NULL, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[2], 
                               Offers_Init = LP_7_Init_extra_2, 
                               Profit_Init = LP_7_Profit_extra_2, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 2, seed = seed, 
                               eps = LP_eps)

## Period 1
# Initial local search - Period 1
ILS_4_gs_extra_1 <- sapply((-Tries):Tries, function(ils)
  ILS_grid_search(ILS_grid = T_Offers, Step_Shift = Step_Shift, ILS_Ind = ils, 
                 Step_Num = Step_Num, sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1,
                 Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB, P_1 = P_Old_Base, 
                 C_1 = Costs_1, Offer_i_t = list(LP_4_extra_3[2, ], LP_7_extra_2[2, ]), 
                 Data_i_t = list(Data_i_3, Data_i_2), GPS_Model = GPS_XGB, 
                 Remain_i_t = 1, t = 1))
Summary_4_gs_extra_1 <- matrix(unlist(ILS_4_gs_extra_1[2, ]), nrow = 2, ncol = 2 * (2 * Tries + 1))
ILS_4_gs_extra_1 <- matrix(unlist(ILS_4_gs_extra_1[1, ]), nrow = N_Indiv_Adj, ncol = 2 * (2 * Tries + 1))
# Closest, optimal (and feasible) initial values - Period 1
Iter_4_gs_extra_1 <- which(Summary_4_gs_extra_1[1, ] <= (alpha_t[1] + NLP_eps))
Ind_4_gs_extra_1 <- Iter_4_gs_extra_1[which.min(Summary_4_gs_extra_1[2, Iter_4_gs_extra_1])]
LP_4_Init_extra_1 <- ILS_4_gs_extra_1[, Ind_4_gs_extra_1]
LP_4_Profit_extra_1 <- Summary_4_gs_extra_1[2, Ind_4_gs_extra_1]
# LP optimizations - Period 1
LP_4_extra_1 <- LP_opt_wrapper(Offer_i_t = list(LP_4_extra_3[2, ], LP_7_extra_2[2, ]), 
                               P_Old_i_1 = P_Old_Base, Data_i_t = list(Data_i_3, Data_i_2), 
                               GPS_Model = GPS_XGB, Remain_i_t = 1, Data_i_1 = Data_i_1, 
                               sd_GPS = sd_GPS_XGB, Mean = Mean_1, Normalizer = Normalizer_1, 
                               Churn_Model = Churn_XGB, RHS_con = LP_RHS, Sign_con = LP_sign, 
                               Sparse_con = LP_con_sparse_L, alpha = alpha_t[1], 
                               Offers_Init = LP_4_Init_extra_1, 
                               Profit_Init = LP_4_Profit_extra_1, Num = Num_L, 
                               Scales = Scales_L, LB = as.numeric(Q_bins[1]), 
                               UB = as.numeric(tail(Q_bins, 1)), t = 1, seed = seed, 
                               eps = LP_eps)

### Additional LP implied churn and profit
# Churn
Churn_LP_Add_1 <- NLP_Churn_i_1(Offers = LP_4_extra_1[2, ], sd_GPS = sd_GPS_XGB, Data_i_1 = Data_i_1, 
                            Mean = Mean_1, Normalizer = Normalizer_1, Churn_Model = Churn_XGB)
Churn_LP_Add_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + LP_4_extra_1[2, ]), 
                            GPS_Model = GPS_XGB, Offers = LP_7_extra_2[2, ], 
                            sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, Churn_Model = Churn_XGB)
Churn_LP_Add_3 <- NLP_Churn_i_t(P_Old = P_Old_Base * (1 + LP_4_extra_1[2, ]) * (1 + LP_7_extra_2[2, ]), 
                            GPS_Model = GPS_XGB, Offers = LP_4_extra_3[2, ], 
                            sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, Churn_Model = Churn_XGB)
# Profit
Profit_LP_Add <- list()
Profit_LP_Add[[1]] <- (1 - Churn_LP_Add_1) * (P_Old_Base * (1 + LP_4_extra_1[2, ]) - Costs_1)
Profit_LP_Add[[2]] <- (1 - Churn_LP_Add_1) * (1 - Churn_LP_Add_2) * 
  (P_Old_Base * (1 + LP_4_extra_1[2, ]) * (1 + LP_7_extra_2[2, ]) - Costs_2)
Profit_LP_Add[[3]] <- (1 - Churn_LP_Add_1) * (1 - Churn_LP_Add_2) * (1 - Churn_LP_Add_3) * 
  (P_Old_Base * (1 + LP_4_extra_1[2, ]) * (1 + LP_7_extra_2[2, ]) * (1 + LP_4_extra_3[2, ]) - Costs_3)

##### Additional NLP multi-period optimization #####
### NLP optimization algorithm
NLP_Opts <- list('algorithm' = 'NLOPT_LN_SBPLX', 'ftol_rel' = 0, 'ftol_abs' = 0,
                 'xtol_rel' = 0, 'xtol_abs' = rep(0, nrow(AMI_i)),
                 'maxeval' = 350000, 'print_level' = 0)

### NLP initial values of previous optimization
NLP_Add_Init <- rep(NA, nrow(Results_Offers_1) * Periods_Max)
NLP_Add_Init[Ind_t_1] <- LP_4_extra_1[2, ]
NLP_Add_Init[Ind_t_2] <- LP_7_extra_2[2, ]
NLP_Add_Init[Ind_t_3] <- LP_4_extra_3[2, ]

### NLP optimization
# Ensure reproducibility by setting the seed for the random number generator
set.seed(seed)
# Retrieve NLP solution
NLP_Add_opt <- nloptr(x0 = NLP_Add_Init, eval_f = NLP_obj_con,
                      lb = NLP_LB_G, ub = NLP_UB_G, opts = NLP_Opts)

### NLP results
NLP_Add_sol <- NLP_Add_opt$solution
NLP_Add_profit <- NLP_Add_opt$objective
NLP_Add_status <- NLP_Add_opt$status

### NLP implied churn and profit
# Churn
Churn_NLP_Add_1 <- NLP_Churn_i_1(Offers = NLP_Add_sol[Ind_t_1], sd_GPS = sd_GPS_XGB, 
                                 Data_i_1 = Data_i_1, Mean = Mean_1, Normalizer = Normalizer_1, 
                                 Churn_Model = Churn_XGB)
Churn_NLP_Add_2 <- NLP_Churn_i_t(P_Old_t = P_Old_Base * (1 + NLP_Add_sol[Ind_t_1]), 
                                 GPS_Model = GPS_XGB, Offers = NLP_Add_sol[Ind_t_2], 
                                 sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_2, 
                                 Churn_Model = Churn_XGB)
Churn_NLP_Add_3 <- NLP_Churn_i_t(P_Old = P_Old_Base * (1 + NLP_Add_sol[Ind_t_1]) * 
                                   (1 + NLP_Add_sol[Ind_t_2]), 
                                 GPS_Model = GPS_XGB, Offers = NLP_Add_sol[Ind_t_3], 
                                 sd_GPS = sd_GPS_XGB, Data_i_t = Data_i_3, 
                                 Churn_Model = Churn_XGB)
# Profit
Profit_NLP_Add <- list()
Profit_NLP_Add[[1]] <- (1 - Churn_NLP_Add_1) * (P_Old_Base * (1 + NLP_Add_sol[Ind_t_1]) - Costs_1)
Profit_NLP_Add[[2]] <- (1 - Churn_NLP_Add_1) * (1 - Churn_NLP_Add_2) * 
  (P_Old_Base * (1 + NLP_Add_sol[Ind_t_1]) * (1 + NLP_Add_sol[Ind_t_2]) - Costs_2)
Profit_NLP_Add[[3]] <- (1 - Churn_NLP_Add_1) * (1 - Churn_NLP_Add_2) * (1 - Churn_NLP_Add_3) * 
  (P_Old_Base * (1 + NLP_Add_sol[Ind_t_1]) * (1 + NLP_Add_sol[Ind_t_2]) * (1 + NLP_Add_sol[Ind_t_3]) - Costs_3)

##### Results #####
### Comparison of optimization offers
cbind(c('LP', 'NLP', 'Additional LP', 'Additional NLP'), 
      c(sum(Profit_LP[[1]]) + sum(Profit_LP[[2]]) + sum(Profit_LP[[3]]),
        sum(Profit_NLP[[1]]) + sum(Profit_NLP[[2]]) + sum(Profit_NLP[[3]]),
        sum(Profit_LP_Add[[1]]) + sum(Profit_LP_Add[[2]]) + sum(Profit_LP_Add[[3]]),
        sum(Profit_NLP_Add[[1]]) + sum(Profit_NLP_Add[[2]]) + sum(Profit_NLP_Add[[3]])))
rbind(c('alpha', 'LP', 'NLP', 'Additional LP', 'Additional NLP'),
      cbind(alpha_t, 
            c(mean(Churn_LP_1), mean(Churn_LP_2), mean(Churn_LP_3)), 
            c(mean(Churn_NLP_1), mean(Churn_NLP_2), mean(Churn_NLP_3)),
            c(mean(Churn_LP_Add_1), mean(Churn_LP_Add_2), mean(Churn_LP_Add_3)),
            c(mean(Churn_NLP_Add_1), mean(Churn_NLP_Add_2), mean(Churn_NLP_Add_3))))

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
Opt_Offer[, 1:Periods_Max] <- cbind(NLP_Add_sol[Ind_t_1], NLP_Add_sol[Ind_t_2], NLP_Add_sol[Ind_t_3])

### Optimal and actual churn
# Pre-allocate output
Opt_Churn <- matrix(rep(NA, N_Indiv_Adj * 2 * Periods_Max), 
                    nrow = N_Indiv_Adj, ncol = 2 * Periods_Max)
# Actual churn
Opt_Churn[, c(Periods_Max + 1:Periods_Max)] <- cbind(Churn_Act_1, Churn_Act_2, Churn_Act_3)
# Optimal churn
Opt_Churn[, 1:Periods_Max] <- cbind(Churn_NLP_Add_1, Churn_NLP_Add_2, Churn_NLP_Add_3)

### Store the final results
