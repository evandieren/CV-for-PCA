# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)


MatrixCompletion <- function(X){
  # args: X matrix containing data
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  # Create bivariate index set for observed data
  n <- nrow(X)
  p <- ncol(X)
  mse5 <- rep(0, p)
  biv_index <- matrix(sample(0:1, size = n*p, replace = T), nrow = n, ncol = p)
  
  # Iterative-hard thresholding algorithm
  tol = 1e-4
  for (r in 1:p) {
    steps <- 0
    M_old <- matrix(rep(1,n*p), ncol=p)
    M_new <- X
    while(norm(M_new-M_old, type = "F")/norm(M_new, type="F") > tol){
      M_old = M_new
      svd_M <- svd(M_old)
      svd_M$d[-(1:r)] <- 0
      M_trunc <- svd_M$u %*% diag(svd_M$d) %*% t(svd_M$v)
      M_new <- X * biv_index + M_trunc * !biv_index
      steps <- steps + 1
    }
    print(paste0("MatrixCompletion: Steps ", steps, ", for rank r ", r))
    mse5[r] <- norm(X-M_new, type = "F")
  }
  return(mse5)
}