# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)


WrongPCA <- function(X, samples){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), samples containing CV-Folds
  # returns: MSE of the CV
  
  K <- dim(samples)[2]
  p <- ncol(X)
  mse <- rep(0, p)
  
  for (r in 1:p) {
    for (k in 1:K) {
      l1 <- length(samples[,k])
      df_k_fold <- X[samples[,k], ]
      df_k <- X[-samples[,k],]
      
      svd_X <- svd(df_k)
      svd_X$d[-(1:r)] <- 0
      X_trunc <- svd_X$u %*% diag(svd_X$d) %*% t(svd_X$v)
      U <- prcomp(X_trunc,rank.=r)$rotation
      
      # Define projection and project data
      P <- U%*%t(U)
      df_k_proj <- P %*% t(df_k_fold)
      
      # Error of estimated and true missing observation
      mse[r] <- sum(sapply(1:l1, function(s){norm(df_k_fold[s,] - df_k_proj[,s], type = "2")^2/l1})) + mse[r]
      
    }
  }
  return(mse)
}