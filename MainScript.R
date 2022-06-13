# --------------------------------------------
# Main Script  
#
# Author: Axel Eurenius Larsson
#
# Title: "On the Value at Risk Forecasting of the Market Risk for 
#         Large Portfolios Based on Dynamic Factor Models
#         with Multivariate GARCH Specifications" 
#
# Project: Masters thesis (Statistics)
# --------------------------------------------

# --------------------------------------------
#### Step 0. Settings
# --------------------------------------------

# Clear global environment
rm(list=ls())

# Set working directory
setwd()

# Set seed
set.seed(980831)

# Output in non-scientific notation
options(scipen = 6, digits = 4)

# Installing Packages
# PANICr (an old package requiring devtools)
require(devtools)
install_version("PANICr", version = "1.0.0", repos = "http://cran.us.r-project.org")

list_of_packages <- c("vars" , # Vector Autoregressions
                      "factoextra", # PCA-estimation
                      "tidyverse", # GGplot2
                      "PANICr", # Deciding number of PC's (ik)
                      "xdcclarge", # cdcc-garch
                      "xts", # time series
                      "Dowd", # Bootstrapping ES
                      "GAS", # Backtesting VaR
                      "ccgarch2", # CC-GARCH MODELS (MIGHT NEED AN OLDER R VERSION FOR THIS)
                      "ggforce", # facet_zoom
                      "psych", # Descriptive statistics
                      "PerformanceAnalytics", # HS-VaR
                      "tikzDevice" # Plotting
)
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, require, character.only = TRUE)

# Reading Functions
source("DFM-functions.R")
source("CCC-functions.R")
source("DCC-functions.R")
source("cDCC-functions.R")
source("Mrisk-functions.R")

# Reading stored Data 
load("Marketcaps.Rda")
load("Valueweights.Rda")
load("dfVAR_CCC.Rda")
load("dfVAR_CCC_sim1.Rda")
load("dfVAR_CCC_sim5.Rda")
load("dfVAR_CCC_sim10.Rda")
load("dfVAR_DCC.Rda")
load("dfVAR_DCC_sim1.Rda")
load("dfVAR_DCC_sim5.Rda")
load("dfVAR_DCC_sim10.Rda")
load("dfVAR_cDCC.Rda")
load("dfVAR_cDCC_sim1.Rda")
load("dfVAR_cDCC_sim5.Rda")
load("dfVAR_cDCC_sim10.Rda")
load("dfPL_ew.Rda")
load("dfPL_vw.Rda")
load("dfPL_sim.Rda")

### Outline
# 0. Settings
# 1. Reading data
# 2. DFM-VaR Algorithm
# 3. Computing P&L vectors for the different portfolios
# 4. Computing the number of breaches in the simulation 

# --------------------------------------------

# --------------------------------------------
# 1. Reading data
# --------------------------------------------

# Reading data from csv-file
Data <- readxl::read_xlsx("Data_master_thesis.xlsx", col_names = T)

# Creating a time series object, removing date column
Data <- xts(Data[,-1], order.by=as.POSIXct(Data$Date))

# Generating log-returns
Data <- apply(log(Data), 2, diff)


# --------------------------------------------
# 2. DFM-VaR-Algorithm
# --------------------------------------------

##### Initiating variables
iw <- 1 # Setting the sample window. Starts at 1 and shifts "right" by 1
iR <- 261 # the length of the sample window,
iN <- ncol(Data) # Number of assets
iT <- nrow(Data) # Number of d?ys in total sample
iF <- iT - iR # Number of days to forecast (first year used as training window)
iM <- 100 # The number of re-samples to compute future scenarios of X*
iB <- 5000 # Number of simulations (the number of random portfolios)
VaR_quantile <- 0.05 # Using 5% for empirical analysis

##### For the simulation study
# Generating a matrix with random portfolio weights
theta_sim <- matrix(NA, nrow = iB, ncol = iN)
set.seed(980831)
for(b in 1:iB){
  wi <- runif(iN,0,1)
  theta_sim[b,]  <- wi/sum(wi)
}

##### Creating storage for VaR forecasts.
# CCC 
CCC_VaR_ew <- c() # Equally weighted portfolio
CCC_VaR_vw <- c()  # Market portfolio
CCC_Var_sim_1 <- matrix(NA, nrow = iB, ncol = iF) # Simulation, a = 0.01
CCC_Var_sim_5 <- matrix(NA, nrow = iB, ncol = iF) # Simulation, a = 0.05
CCC_Var_sim_10 <- matrix(NA, nrow = iB, ncol = iF) # Simulation, a = 0.10
# DCC 
DCC_VaR_ew <- c()
DCC_VaR_vw <- c()
DCC_Var_sim_1 <- matrix(NA, nrow = iB, ncol = iF)
DCC_Var_sim_5 <- matrix(NA, nrow = iB, ncol = iF)
DCC_Var_sim_10 <- matrix(NA, nrow = iB, ncol = iF)
# cDCC 
cDCC_VaR_ew <- c()
cDCC_VaR_vw <- c()
cDCC_Var_sim_1 <- matrix(NA, nrow = iB, ncol = iF)
cDCC_Var_sim_5 <- matrix(NA, nrow = iB, ncol = iF)
cDCC_Var_sim_10 <- matrix(NA, nrow = iB, ncol = iF)

# Creating a dummy for the while-loop
reached_end_of_sample <- FALSE
year <- 2018

# Making one step ahead forecasts using a 250 day rolling window
while (reached_end_of_sample == FALSE) {
  
  # --------------------------------------------
  # STEP 2.A: Variables, initiating data and so on
  # --------------------------------------------
  
  # Initiating the sample window of length iR.
  # This is the number of trading days in 2017 avaliable
  it <- iw+iR # timepoint t (last in sample)
  sample_window <- as.xts(Data[iw:it, ])
  
  # --------------------------------------------
  # STEP 2.B: Estimating Principal Components
  # --------------------------------------------
  
  # Deciding the number of principal components according to Bai and NG (2002)
  ir <- getnfac(x = sample_window, kmax = 10, criteria = "IC1")$ic
  if(ir == 1){
    ir <- 2 # At least 2 PC's
  }
  
  # Estimating principal components
  PC <- fnPC(y = sample_window, nfac = ir)
  
  # --------------------------------------------
  # STEP 2.C: Estimating Lambda, F and epsilon
  # --------------------------------------------
  
  # Lambda (is N by r)
  L_hat <- PC$lambda
  
  # Static factors (is r by t)
  F_hat <- PC$fhat
  colnames(F_hat) <- sapply(seq_along(1:ir), function(x) paste0("PC",x))
  
  # Epsilon at t (N by T)
  e_hat <- PC$ehat
  
  # --------------------------------------------
  # STEP 2.D: Using the estimated static factors, F,
  #         run the VAR. Extract A, Hu_t, u_t
  # --------------------------------------------
  # Running the VAR
  var_est <- VAR(F_hat, p = 1) # VAR is fixed to 1 lag
  
  # Extracting the coefficient estimates A
  A_hat <- as.matrix(unlist(Acoef(var_est)[[1]]))
  
  # Extracting the residuals, Hu_t
  Hu_hat <- resid(var_est)
  
  # Computing H_hat
  H_hat <- prcomp(Hu_hat, center = F, scale. = F)$rotation
  
  # Computing u_hat "backwards" (k x T)
  u_hat <- Hu_hat%*%H_hat
  
  
  # --------------------------------------------
  # STEP 2.E: CCC-GARCH estimation & forecasting, Only using CCC(1,1)-GARCH(1,1)
  # Fixing starting values in optimization so it starts from the same matrix each time, used in CCC, DCC, cDCC
  start_val <- fnstart_vals(ir)
  # --------------------------------------------
  
  # Estimating the CCC(1,1)-GARCH(1,1)-Normal
  set.seed(980831)
  CCC <- estimateCCC(data = u_hat, iniA = start_val[[1]], iniB = start_val[[2]])
  
  # Extracting estimates of z, standardized residuals
  CCC_z <- as.matrix(CCC$z)
  
  # Extracting estimates of omega, alpha and beta
  CCC_params <- matrix(c(CCC$estimates$a,
                         diag(CCC$estimates$A),
                         diag(CCC$estimates$B)),
                       nrow = ir)
  
  # Forecasting D
  CCC_D_forecast <- fnD_fcast(param_mat = CCC_params,
                              epsilon_t = u_hat[iR, ],
                              h_t = CCC$h[iR, ])
  
  # Forecasting R not neccesary, just use estimated R since it is constant.
  CCC_R_forecast <- CCC$estimates$R
  
  # Forecasting H
  CCC_H_forecast <- fnH_fcast(Rforecast = CCC_R_forecast, 
                              Dforecast = CCC_D_forecast)
  
  # --------------------------------------------
  # STEP 2.F: DCC-GARCH estimation & forecasting, Only using DCC(1,1)-GARCH(1,1)
  # --------------------------------------------
  
  # Estimating the DCC(1,1)-GARCH(1,1)-Normal
  set.seed(980831)
  DCC <- estimateDCC(data = u_hat, iniA = start_val[[1]], iniB = start_val[[2]])
  
  # Extracting estimates of z, standardized residuals
  DCC_z <- as.matrix(DCC$z)
  
  # Extracting estimates of omega, alpha and beta (not extracting mu which are the ir first estimates)
  DCC_params <- matrix(DCC$f.stage$pars[(ir+1):(4*ir)], nrow = ir)
  
  # Forecasting D
  DCC_D_forecast <- fnD_fcast(param_mat = DCC_params,
                              epsilon_t = u_hat[iR, ],
                              h_t = DCC$h[iR, ])
  
  # Forecasting R
  DCC_R_forecast <- fnR_fcast_DCC(DCC_a = DCC$s.stage$pars[1],
                                  DCC_b = DCC$s.stage$pars[2],
                                  DCC_e = DCC_z)
  
  # Forecasting H
  DCC_H_forecast <- fnH_fcast(Rforecast = DCC_R_forecast, 
                              Dforecast = DCC_D_forecast)
  
  # --------------------------------------------
  # STEP 2.G: cDCC-GARCH-estimation & forecasting
  # --------------------------------------------
  
  # Estimating cDCC model
  set.seed(980831)
  cDCC <- estimateCDCC(data = u_hat, iniA = start_val[[1]], iniB = start_val[[2]])
  
  # Extracting estimates of z, standardized residuals
  cDCC_z <- as.matrix(cDCC$z)
  
  # Extracting estimates of omega, alpha and beta (not extracting mu which are the ir first estimates)
  cDCC_params <- matrix(cDCC$f.stage$pars[(ir+1):(4*ir)], nrow = ir)
  
  # Forecasting the D matrix
  cDCC_D_forecast <- fnD_fcast(param_mat = cDCC_params,
                               epsilon_t = u_hat[iR, ],
                               h_t = cDCC$h[iR, ])
  
  # Forecasting the R matrix
  cDCC_R_forecast <- fnR_fcast_cDCC(cDCC_a = cDCC$s.stage$pars[1],
                                    cDCC_b = cDCC$s.stage$pars[2],
                                    cDCC_e = cDCC_z)
  
  # Forecasting the H matrix
  cDCC_H_forecast <- fnH_fcast(Rforecast = cDCC_R_forecast, 
                               Dforecast = cDCC_D_forecast)
  
  # --------------------------------------------
  # STEP 2.H: Building scenarios for P&L and
  #         Estimating DFM-VaR and Bootstrapped-ES
  # --------------------------------------------
  
  ### CCC-estimations, Normal distribution
  CCC_risk <- risk_estimation(sample_length = iR,
                              weight = portfolio_weights,
                              H_forecast = CCC_H_forecast,
                              z_forecast = CCC_z,
                              Quantile = VaR_quantile,
                              year = year)
  
  CCC_VaR_ew[iw] <- CCC_risk$VaR_ew
  CCC_VaR_vw[iw] <- CCC_risk$VaR_vw
  CCC_Var_sim_1[,iw] <- CCC_risk$Var_sim1
  CCC_Var_sim_5[,iw] <- CCC_risk$Var_sim5
  CCC_Var_sim_10[,iw] <- CCC_risk$Var_sim10
  
  
  ### DCC-estimations, Normal distribution
  DCC_risk <- risk_estimation(sample_length = iR,
                              weight = portfolio_weights,
                              H_forecast = DCC_H_forecast,
                              z_forecast = DCC_z,
                              Quantile = VaR_quantile,
                              year = year)
  
  DCC_VaR_ew[iw] <- DCC_risk$VaR_ew
  DCC_VaR_vw[iw] <- DCC_risk$VaR_vw
  DCC_Var_sim_1[,iw] <- DCC_risk$Var_sim1
  DCC_Var_sim_5[,iw] <- DCC_risk$Var_sim5
  DCC_Var_sim_10[,iw] <- DCC_risk$Var_sim10
  
  
  
  ### cDCC-estimations
  cDCC_risk <- risk_estimation(sample_length = iR,
                               weight = portfolio_weights,
                               H_forecast = cDCC_H_forecast,
                               z_forecast = cDCC_z,
                               Quantile = VaR_quantile,
                               year = year)
  
  cDCC_VaR_ew[iw] <- cDCC_risk$VaR_ew
  cDCC_VaR_vw[iw] <- cDCC_risk$VaR_vw
  cDCC_Var_sim_1[,iw] <- cDCC_risk$Var_sim1
  cDCC_Var_sim_5[,iw] <- cDCC_risk$Var_sim5
  cDCC_Var_sim_10[,iw] <- cDCC_risk$Var_sim10
  
  # --------------------------------------------
  # STEP 2.I: If we have reached the end of the sample:
  #        End the loop. Otherwise, move one the
  #        sample window 1 step forward in time.
  # --------------------------------------------
  # Changing year index for Value-Weighted portfolio
  if(as.numeric(substr(rownames(Data)[it],0,4)) != year){
    year <- year + 1
  }
  
  if(it == nrow(Data)){
    reached_end_of_sample = TRUE
  } else {
    iw <- iw+1
  }
}

# --------------------------------------------
# 3. Computing P&L vectors for the different portfolios
# --------------------------------------------

# Equally weighted portfolio
theta_equalweight <- matrix(rep(1/iN, iN), nrow = 1, ncol = iN)
realized_PL_ew <- Data[(iR+1):(iR+iF),]%*%t(theta_equalweight)

# Value weighted portfolio
realized_PL_vw2018 <- Data[rownames(Data)[as.numeric(substr(rownames(Data),0,4)) == 2018],]%*%t(as.matrix(Valueweights[1,]))
realized_PL_vw2019 <- Data[rownames(Data)[as.numeric(substr(rownames(Data),0,4)) == 2019],]%*%t(as.matrix(Valueweights[2,]))
realized_PL_vw2020 <- Data[rownames(Data)[as.numeric(substr(rownames(Data),0,4)) == 2020],]%*%t(as.matrix(Valueweights[3,]))
realized_PL_vw2021 <- Data[rownames(Data)[as.numeric(substr(rownames(Data),0,4)) == 2021],]%*%t(as.matrix(Valueweights[4,]))
realized_PL_vw <- matrix(c(realized_PL_vw2018[2:length(realized_PL_vw2018)],realized_PL_vw2019,realized_PL_vw2020,realized_PL_vw2021),
                         nrow = iF, ncol = 1)
rownames(realized_PL_vw) <- rownames(realized_PL_ew)

# Simulated portfolio
realized_PL_sim <- theta_sim%*%t(Data[(iR+1):nrow(Data), ])

# --------------------------------------------
# 4. Computing the number of breaches in the simulation 
# --------------------------------------------
breach_CCC_1 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_CCC_sim1)
breach_CCC_5 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_CCC_sim5)
breach_CCC_10 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_CCC_sim10)
breach_DCC_1 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_DCC_sim1)
breach_DCC_5 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_DCC_sim5)
breach_DCC_10 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_DCC_sim10)
breach_cDCC_1 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_cDCC_sim1)
breach_cDCC_5 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_cDCC_sim5)
breach_cDCC_10 <- fnVar_breaches_sim(realized_PL_sim, dfVAR_cDCC_sim10)
