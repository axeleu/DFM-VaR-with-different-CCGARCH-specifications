### Estimating DFM-VaR
risk_estimation <- function(sample_length, weight, H_forecast, z_forecast, Quantile, year){
  
  # Portfolio weights, equally weighted
  theta_ew <- matrix(rep(1/iN, iN), nrow = 1, ncol = iN)
  
  # Portfolio weights, market portfolio, weights changed yearly
  if(year == 2018) { 
    theta_vw <- as.matrix(Valueweights[1,])
  }else if(year == 2019) {
    theta_vw <- as.matrix(Valueweights[2,])
  }else if(year == 2020) {
    theta_vw <- as.matrix(Valueweights[3,])
  }else if(year == 2021) {
    theta_vw <- as.matrix(Valueweights[4,])
  }
  
  # For storage, profit & loss vector
  PL_ew <- c() # Equally weighted
  PL_vw <- c() # Value Weighted 
  PL_sim <- matrix(NA, nrow = iB, ncol = iM) # Randomized 
  

  #-----------------------------------------
  # Value-At-Risk Estimation
  #-----------------------------------------
  
  # Square root of matrix such that A = BB
  # Using eigenvalue decomposition (requires PSD matrix)
  eigen_H <- eigen(H_forecast)
  sqrt_H <- eigen_H$vectors %*% diag(sqrt(eigen_H$values)) %*% solve(eigen_H$vectors)
  
  # Random sample
  set.seed(980831)
  sample_w <- sample(1:sample_length, iM)
  
  for (iter in seq_along(sample_w)) {
    
    # Building scenarios for X^* at t+1
    X_t1 <- L_hat%*%(A_hat%*%F_hat[sample_length+1,] + H_hat%*%sqrt_H%*%z_forecast[sample_w[iter], ]) + t(e_hat[sample_w[iter], ])
    
    # Scenarios of Profit and Loss (P&L) at t+1 using equally weighted portfolio
    PL_ew[iter] <- theta_ew%*%X_t1
    
    # Scenarios of Profit and Loss (P&L) at t+1 using value weighted portfolio
    PL_vw[iter] <- theta_vw%*%X_t1
    
    # theta_sim is of dim 5000*132
    # X_t1 is of dim 132*(100)
    # PL_sim is of dim 5000*100 and consists of 
    # the P&L for each randomly weighted portfolio 
    # and each simulated X*.
    PL_sim[,iter] <- theta_sim%*%X_t1
  }

  # Computing VaR for the equally weigted portfolio
  VaR_ew <- quantile(PL_ew, Quantile)
  
  # Computing VaR for the market portfolio
  VaR_vw <- quantile(PL_vw, Quantile)
  
  # Computing VaR for the simulated random portfolios
  # Using 0.01, 0.05, 0.10 as confidence level
  Var_sim1  <- apply(PL_sim, 1, function(x) quantile(x, 0.01))
  Var_sim5  <- apply(PL_sim, 1, function(x) quantile(x, 0.05))
  Var_sim10 <- apply(PL_sim, 1, function(x) quantile(x, 0.10))
  
  # Returning
  return(list("VaR_ew" = VaR_ew,
              "VaR_vw" = VaR_vw,
              "Var_sim1" = Var_sim1,
              "Var_sim5" = Var_sim5,
              "Var_sim10" = Var_sim10))
}


# Empirical application breaches function (one portfolio)
fnVar_breaches_emp<- function(PL, VaR) {
  no_breaches <- sum(PL < VaR)
  breach_indicator <- ifelse(PL < VaR, 1, 0)
  return(list("No.Breaches" = no_breaches,
         "Breach_ind" = breach_indicator))
}


# Simulation study breaches function (B portfolios)
fnVar_breaches_sim <- function(realized_PL_sim, simulated_VaR) {
  no_breaches <- c()
  for (b in 1:iB) {
    no_breaches[b] <- sum(realized_PL_sim[b, ] < simulated_VaR[b, ])
  }
  mean_breaches <- mean(no_breaches)
  var_breaches <- var(no_breaches)

  
  return(list("no_breaches" = no_breaches,
              "mean_breach" = mean_breaches,
              "var_breach" = var_breaches))
}