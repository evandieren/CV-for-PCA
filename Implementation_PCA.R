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


WrongPCAImproved <- function(X, samples){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), samples containing CV-Folds
  # returns: MSE of the CV
  
  K <- dim(samples)[2]
  p <- ncol(X)
  mse1 <- rep(0, p)
  
  # Split data into missing and observed data
  split <- sample(1:p, floor(p/2))
  
  for (r in 1:p) {
    for (k in 1:K) {
      l1 <- length(samples[,k]) # number of elements in the fold
      
      df_k <- X[-samples[,k],] ## all observations except the fold
      mu <- colMeans(df_k) # mu without fold
      svd_sigma <- svd(cov(df_k)) # eigen here
      svd_sigma$d[-(1:r)] <- 0
      sigma_trunc <- svd_sigma$u %*% diag(svd_sigma$d) %*% t(svd_sigma$v) # truncating the covariance matrix
      
      df_k_fold <- X[samples[,k],]
      df_k_fold_miss <- as.matrix(df_k_fold[,split])
      df_k_fold_obs <- as.matrix(df_k_fold[,-split])
      
      # Estimate est_x_miss for df_k_fold_miss
      mu_miss <- mu[split]
      mu_obs <- mu[-split]
      sigma_miss_obs <- sigma_trunc[split, -split]
      sigma_obs_obs <- sigma_trunc[-split, -split]
      
      est_x_miss <- sapply(1:l1, function(n){mu_miss + sigma_miss_obs %*% ginv(sigma_obs_obs, tol = 1e-20) %*% (df_k_fold_obs[n,]-mu_obs)})
      #print(est_x_miss)
      #mse1[r] <- sum(sapply(1:l1, function(s){norm(est_x_miss[[s]]-df_k_fold_miss[s,], type = "2")^2/l1})) + mse1[r]
      #mse1[r] <- sum(sapply(1:l1, function(s){sum((est_x_miss[[s]]-df_k_fold_miss[s,])^2)}))/l1 + mse1[r]
      mse1[r]  <-  sum((est_x_miss - t(df_k_fold_miss))^2) + mse1[r]
      
    }
  }
  return(mse1)
}

KDEApproach <- function(X){
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
    M_old <- matrix(rep(1, size = n*p), nrow = n, ncol = p)
    M_new <- X
    while(norm(M_new-M_old, type = "F")>tol){ # add /norm(M_new,type="F") ?
      M_old = M_new
      svd_M <- svd(M_old)
      svd_M$d[-(1:r)] <- 0
      M_trunc <- svd_M$u %*% diag(svd_M$d) %*% t(svd_M$v)
      M_new <- X * biv_index + M_trunc * !biv_index
      steps <- steps + 1
    }
    print(paste0("MatrixCompletion: Steps ", steps, ", for rank r ", r))
    mse5[r] <- norm(X-M_new, type = "F")^2
  }
  return(mse5)
}

set.seed(11)
n <- 100
p <- 10
X <- array(rnorm(121,100,1),c(n,p))
svd_X <- svd(X)
svd_X$d[-(1:p)] <- 0
X_tronc <- svd_X$u %*% diag(svd_X$d) %*% t(svd_X$v)
samples <- matrix(sample(1:n),ncol=5)

#out <- MatrixCompletion(X)
#out <- KDEApproach(X)
#out <- WrongPCAImproved(X_tronc,samples)

# Testing Algorithms - seem to work so far

# n <- 100
# p <- 10
# K <- 10
# mean_df <- rep(0, p)
# mat <- matrix(rnorm(p*p, mean = 0, sd = 12), nrow = p, ncol = p)
# df <- rmvnorm(n = n, mean = mean_df, sigma = mat %*% t(mat))
# samples <- matrix(sample(1:n),ncol=K)
# 
# WrongPCA(df, samples)
# WrongPCAImproved(df, samples)
# MatrixCompletion(df)
# KDEApproach(df)










# Questions:

# 25.11.22
# How is r choosen? -> PCA first
# How missing and observed data split? -> randomly, but rule of thumb is half half
# MSE on x_miss? -> Yes
# Less informations in Method 1? -> No, different computations of mu and sigma
# Method 1: Truncating Missed and Observed data accordingly? -> No
# Method 2: Still p dimensions? -> Yes
# Method 3: How to cross-validate? Different error measure of covariance similiarity








