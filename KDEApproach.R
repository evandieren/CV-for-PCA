# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)


KDEApproach <- function(X, samples = NaN){
  # args: X matrix containing data
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  p <- ncol(X)
  n <- nrow(X)
  mse4 <- rep(0,p)
  
  for (r in 1:p) {
    eigen_sigma <- eigen(cov(X))
    eigen_sigma$values[-(1:r)] <- 0
    sigma_trunc <- eigen_sigma$vectors %*% diag(eigen_sigma$values) %*% t(eigen_sigma$vectors)
    
    mse4[r] <- mse4[r] + norm(sigma_trunc,type="F")^2
    
    for (i in 1:n){
      eigen_sigma_n <- eigen(cov(X[-i,]))
      eigen_sigma_n$values[-(1:r)] <- 0
      sigma_trunc_n <- eigen_sigma_n$vectors %*% diag(eigen_sigma_n$values) %*% t(eigen_sigma_n$vectors)
      mse4[r] <- mse4[r] - 2/n*sum(sigma_trunc_n*(X[i,]%*%t(X[i,])))
    }
    
  }
  return(mse4)
}

#set.seed(11)
n <- 100
p <- 5
r <- 3
X <- array(rnorm(n*p,0,1),c(n,p))
svd_X <- svd(X)
svd_X$d[-(1:r)] <- 0
X_tronc <- svd_X$u %*% diag(svd_X$d) %*% t(svd_X$v)
X_tronc_noise <- X_tronc + rmvnorm(100, mean=rep(0, p), sigma=svd_X$d[r]*0.1*diag(p))
#samples <- matrix(sample(1:n),ncol=5)


out_tronc <- KDEApproach(X_tronc)
out_tronc_noise <- KDEApproach(X_tronc_noise)
par(mfrow=c(1,2))
plot(1:p, out_tronc, "l", col=1)
plot(1:p, out_tronc_noise, "l", col=2)
