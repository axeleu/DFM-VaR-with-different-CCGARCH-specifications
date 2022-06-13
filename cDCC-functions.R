fnR_fcast_cDCC <- function(cDCC_a, cDCC_b, cDCC_e){
  # Following Aielli (2013)
  # Estimated parameters from cDCC estimation from ccgarcg2-package
  alpha <- cDCC_a
  beta <- cDCC_b
  epstilde <- cDCC_e
  
  #----------------------------------------------------------------------
  ##### Definition 3.3 (cDCC conditional estimator of S) in Aielli (2013)
  #----------------------------------------------------------------------
  
  # Creating initial values
  initial_value <- c()
  for(iter in 1:ir){
    initial_value[iter] <- (1 - alpha - beta)
  }

  # Start by estimating sqrt of Q tilde star (Qtildestarsqrt)
  Qtildestarsqrt <- list()
  Qtildestarsqrt[[1]] <- diag(sqrt(initial_value), ir, ir)
  # Elementwise estimation, until iR+1 which is the forecast
  for (k in 2:iR) {
    qtilde <- c()
    for (asset in 1:ir) {
      # Equation 18 in Aielli
      qtilde[asset] <- (1-alpha-beta) + alpha*(epstilde[k-1,asset]**2)*Qtildestarsqrt[[k-1]][asset,asset] + beta*Qtildestarsqrt[[k-1]][asset,asset] 
    }
    Qtildestarsqrt[[k]] <- diag(sqrt(qtilde), ir, ir)
  }

  # Computing Q%*%eps%*%t(eps)%*%t(Q)
  Qtildesqrt_epstilde <- lapply(1:iR, 
                                function(i,Q,e){Q[[i]]%*%e[i,]%*%t(e[i,])%*%t(Q[[i]])}, # Second Q transpose to guarantee PD. 
                                Q=Qtildestarsqrt,e=epstilde) 
  
  # S_hat = (1/T)*sum(Q_t%*%eps_t%*%t(eps_t)%*%Q_t)
  S_hat <- Reduce("+", Qtildesqrt_epstilde) / length(Qtildesqrt_epstilde)
  
  #----------------------------------------------------------------------
  #### Equation 11 in Aielli (2013)
  #----------------------------------------------------------------------
  Q <- list()
  Q[[1]] <- (1-alpha-beta)*S_hat + beta*S_hat # Initial values for Q ([1-a-b]S + bS = S - aS - bS + bS = (1-a)S)
  
  Qstarsqrt <- list()
  Qstarsqrt[[1]] <- diag(diag(Q[[1]])) # Initial values for Qstarsqrt
  
  # Estimating Q until iR+1 which is the forecast
  for (t in 2:(iR+1)) {
    Q[[t]] <-  (1 - alpha - beta)*S_hat + alpha*(Qstarsqrt[[t-1]]%*%epstilde[t-1,]%*%t(epstilde[t-1,])%*%t(Qstarsqrt[[t-1]])) + beta*Q[[t-1]]
    Qstarsqrt[[t]] <- diag(sqrt(diag(Q[[t]])))
    
  }
  
  # Forecasting R using forecasted Q.
  mR_cDCC_forecast <- diag(1/(sqrt(diag(Q[[iR+1]]))))%*%Q[[iR+1]]%*%diag(1/(sqrt(diag(Q[[iR+1]]))))
  
  return(mR_cDCC_forecast)
}

