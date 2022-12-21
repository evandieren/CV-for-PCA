# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)


MatrixCompletion <- function(X, samples = NaN){
  # args: X matrix containing data
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  # Create bivariate index set for observed data
  n <- nrow(X)
  p <- ncol(X)
  mse5 <- rep(0, p)
  
  percent_miss <- floor(p/2)

  matrix_miss <- matrix(data=0,n,percent_miss)
  for (j in 1:n)
    matrix_miss[j,] <- sample(1:p,percent_miss)
  
  miss_X <- X
  for (i in 1:n){
    miss_X[i,matrix_miss[i,]] <- NA
  }
  
  mu_obs <- colMeans(miss_X,na.rm = TRUE)
  
  # Iterative-hard thresholding algorithm
  for (r in 1:p) {
    X0 <- matrix(rnorm(n=n*p, mean=mu_obs, sd=1), ncol=p)
    #X0 <- rmvnorm(n,mu_obs,cov(X))
    
    M_new <- MatrixCompletion_only(miss_X,X0,r)
    
    for (i in 1:n){
      mse5[r] <- mse5[r] + sum((X[i,matrix_miss[i,]]-M_new[i,matrix_miss[i,]])^2)
    }
    
    mse5[r] <- mse5[r]/(n*percent_miss)
  }
  return(mse5)
}

# ALMOND: extract only core function for clearity --------------------------

MatrixCompletion_only <- function(X, X0, r, tol = 1e-4, maxsteps = 1000){
  # args: X: matrix containing data - missing values are NA
  # X0: complete matrix of the same format as X, initial value for algorithm
  # tol, maxsteps: error tolerance and maximum steps, so algo doesn't run forever
  # returns: Completed matrix
  
  # Create bivariate index set for observed data
  n <- nrow(X)
  p <- ncol(X)
  biv_index <- is.na(X)
  
  steps <- 1
  crit <- TRUE
  M_new <- X0
  
  # Iterative-hard thresholding algorithm
  while(crit > tol){
      M_old <- M_new
      svd_M <- svd(M_old)
      svd_M$d[-(1:r)] <- 0
      M_trunc <- svd_M$u %*% diag(svd_M$d) %*% t(svd_M$v)
      M_new <- ifelse(is.na(X), M_trunc, X)
      steps <- steps + 1
      crit <- norm(M_new-M_old, type="F")/norm(M_new, type="F")
      if(steps == maxsteps) {
        warning(paste("Stopped after iteration", steps, 
                      "with stopping criterion", crit, "> tol."))
        break
      }
    }
    print(paste0("MatrixCompletion: Steps ", steps, ", for rank r ", r))
  return(M_new)
}


