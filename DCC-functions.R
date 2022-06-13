### Forecasting R
fnR_fcast_DCC <- function(DCC_a, DCC_b, DCC_e) {

  # Estimated parameters from cDCC estimation from ccgarcg2-package
  alpha <- DCC_a
  beta <- DCC_b
  epstilde <- DCC_e
  
  # Computing S, second moment of standardized resiudals
  S <- cov(epstilde)
  
  # Computing Q, loop goes until iR+1 which is the forecasted value
  Q <- list()
  Q[[1]] <- (1-alpha-beta)*S + beta*S
  for (it in 2:(iR+1)) {
    Q[[it]]<- (1-alpha-beta)*S + (alpha*epstilde[it-1,]%*%t(epstilde[it-1,])) + beta*Q[[it-1]]
  }
  
  
  mR_forecast_DCC <- diag(1/(sqrt(diag(Q[[iR+1]]))))%*%Q[[iR+1]]%*%diag(1/(sqrt(diag(Q[[iR+1]]))))
  return(mR_forecast_DCC)
}


