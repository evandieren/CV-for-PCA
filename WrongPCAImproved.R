# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)


WrongPCAImproved <- function(X, samples){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), samples containing CV-Folds
  # returns: MSE of the CV
  
  K <- dim(samples)[2]
  p <- ncol(X)
  mse1 <- rep(0, p)
  
  # Split data into missing and observed data
  split <- sample(1:p, floor(p/2)) # 2 - 4
  
  for (r in 1:p) {
    for (k in 1:K) {
      l1 <- length(samples[,k]) # number of elements in the fold
      
      df_k <- X[-samples[,k],] ## all observations except the fold
      mu <- colMeans(df_k) # mu without fold
      eigen_sigma <- eigen(cov(df_k)) # eigen here
      sum_eigen <- sum(eigen_sigma$d[-(1:r)])
      eigen_sigma$d[-(1:r)] <- 0
      sigma_trunc <- eigen_sigma$u %*% diag(eigen_sigma$d) %*% t(eigen_sigma$v) + sum_eigen/p*diag(p)# truncating the covariance matrix
      
      df_k_fold <- X[samples[,k],]
      df_k_fold_miss <- as.matrix(df_k_fold[,split])
      df_k_fold_obs <- as.matrix(df_k_fold[,-split])
      
      # Estimate est_x_miss for df_k_fold_miss
      mu_miss <- mu[split]
      mu_obs <- mu[-split]
      sigma_miss_obs <- sigma_trunc[split, -split]
      sigma_obs_obs <- sigma_trunc[-split, -split]
      
      est_x_miss <- sapply(1:l1, function(n){mu_miss + sigma_miss_obs %*% ginv(sigma_obs_obs) %*% (df_k_fold_obs[n,]-mu_obs)})
      #print(est_x_miss)
      #mse1[r] <- sum(sapply(1:l1, function(s){norm(est_x_miss[[s]]-df_k_fold_miss[s,], type = "2")^2/l1})) + mse1[r]
      #mse1[r] <- sum(sapply(1:l1, function(s){sum((est_x_miss[[s]]-df_k_fold_miss[s,])^2)}))/l1 + mse1[r]
      mse1[r]  <-  sum((est_x_miss - t(df_k_fold_miss))^2) + mse1[r]
      
    }
  }
  return(mse1)
}







# Questions:

# 25.11.22
# How is r choosen? -> PCA first
# How missing and observed data split? -> randomly, but rule of thumb is half half
# MSE on x_miss? -> Yes
# Less informations in Method 1? -> No, different computations of mu and sigma
# Method 1: Truncating Missed and Observed data accordingly? -> No
# Method 2: Still p dimensions? -> Yes
# Method 3: How to cross-validate? Different error measure of covariance similiarity








