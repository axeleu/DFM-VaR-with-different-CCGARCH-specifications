# Estimating PC
fnPC <- function(y, nfac) {
  
  bigt <- dim(y)[1]
  
  bign <- dim(y)[2]
  
  covmat <- crossprod(y) / (bigt-1)
  
  eig <- eigen(covmat)
  
  e_vec <- eig$vectors
  
  lambda <- sqrt(bign)*e_vec[,1:nfac]
  
  fhat <- y %*% (lambda/bign)
  
  ehat <- y - tcrossprod(fhat, lambda)
  
  return(list(ehat = ehat, fhat = fhat, lambda = lambda))
}