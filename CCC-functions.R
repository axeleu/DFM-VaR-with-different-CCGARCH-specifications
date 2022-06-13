#### CCC-Estimation
# Step 1. Estimate CCC, Obtain the correlation matrix R
# Step 2. Forecast the diagonal matrix D with h_it on the diagonal
# Step 3. Since Ht = Dt R Dt => Ht+1 = Dt+1 R Dt+1

### Forecasting D matrix (1 day ahead)
fnD_fcast <- function(param_mat, epsilon_t, h_t){
  # param_mat is ordered as: "omega", "alpha", "beta" by columns and PC by rows. 
  omega <- param_mat[,1]
  alpha <- param_mat[,2]
  beta <- param_mat[,3]
  
  h_forecast <- c()
  for (asset in 1:ir) {
    h_forecast[asset] <- omega[asset] + alpha[asset]*(epsilon_t[asset]**2) + beta[asset]*matrix(h_t)[asset]
  }
  
  # Storing D matrix as sqrt of sigma, diagonal elements
  D_forecast <- diag(sqrt(h_forecast))
  
  # Returning D
  return(D_forecast)
}


# Forecasting the H matrix
fnH_fcast <- function(Rforecast, Dforecast){
  return(Dforecast%*%Rforecast%*%Dforecast)
}

# Fixing starting values for estimation. 
fnstart_vals <- function(ir) {
  set.seed(980831)
  ndim <- c(ir)
  ub <- runif(ndim, min=0.0001, max=0.04)
  iniA <- matrix(runif(ndim^2, min=0, max=ub[sample(1:ndim, 1)]), ndim, ndim)
  iniB <- matrix(runif(ndim^2, min=-0.004, max=ub[sample(1:ndim, 1)]), ndim, ndim)
  diag(iniA) <- round(runif(ndim, min=0.04, max=0.05), 4)
  diag(iniB) <- round(runif(ndim, min=0.8, max=0.9), 4)
  return(list(iniA, iniB))
}
