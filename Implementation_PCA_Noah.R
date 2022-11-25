# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)


df_df <- data.frame("obs1"=c(5,43,2,3,1,4,8,4,5,9), "obs2"=c(5,43,2,3,1,4,8,4,5,9), "obs3"=c(1,2,5,3,4,5,6,1,5,7), "obs4"=c(1,7,6,3,4,5,6,1,7,10))
df <- as.matrix(df_df)


# Create Multivariate Random Variables
sigma <- matrix(c(5,43,2,3,1,4,8,4,5,9,5,43,2,3,1,4,8,4,5,9,5,43,2,3,1), nrow = 5, ncol = 5)
df1 <- rmvnorm(n = 20, mean = rep(0, nrow(sigma)), sigma = sigma %*% t(sigma))


# R built-in function for PCA
pca3 <- prcomp(df1, rank. = 5)
summary(pca3) # overview of the extent a variable explains the data
pca3$sdev # array sorted by Principal components that explain % of standard deviation


WrongPCA <- function(X, pca_vec, K){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), K amount of CV-Folds
  # returns: MSE of the CV
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  index <- c(1:n)
  K <- K
  samples <- vector("list", K)
  
  for (i in 1:K) {
    if (i < K) {
      fold <- sample(index, n/K)
      samples[[i]] <- fold
      index <- index[-which(index %in% fold)]
    }
    else{
      samples[[i]] <- index
    }
  }
  
  mse <- rep(0, K)
  
  for (k in 1:K) {
    l1 <- length(samples[[k]])
    df_k <- X[samples[[k]], ]
    
    svd_X <- svd(X)
    X_trunc <- svd_X$u[,pca_vec] %*% diag(svd_X$d[pca_vec]) %*% t(svd_X$v)[pca_vec,]

    # Define projection and project data
    P <- X_trunc %*% solve(t(X_trunc) %*% X_trunc, tol = 1e-20) %*% t(X_trunc)
    X_proj <- P %*% X
    df_k_proj <- X_proj[samples[[k]],]
    
    # Error of estimated and true missing observation
    for(s in 1:l1){
      mse[k] <- mse[k] + norm(df_k[s,] - df_k_proj[s,], type = "2")^2/l1/K
    }
  }
  
  return(sum(mse))

}

WrongPCA(df1, 1:3, 5)


WrongPCAImproved <- function(X, pca_vec, K){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), K amount of CV-Folds
  # returns: MSE of the CV
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  index <- c(1:n)
  K <- K
  samples <- vector("list", K)
  
  for (i in 1:K) {
    if (i < K) {
      fold <- sample(index, n/K)
      samples[[i]] <- fold
      index <- index[-which(index %in% fold)]
    }
    else{
      samples[[i]] <- index
    }
  }
  
  # Split data into missing and observed data
  split <- sample(1:p, floor(p/2))
  
  mse1 <- rep(0, K)
  
  for (k in 1:K) {
    l1 <- length(samples[[k]])
    
    df_k <- X[-samples[[k]],]
    mu <- colMeans(df_k)
    svd_sigma <- svd(cov(df_k))
    sigma_trunc <- svd_sigma$u[,pca_vec] %*% diag(svd_sigma$d[pca_vec]) %*% t(svd_sigma$v)[pca_vec,]
    
    df_k_fold <- X[samples[[k]],]
    df_k_fold_miss <- as.matrix(df_k_fold[,split])
    df_k_fold_obs <- as.matrix(df_k_fold[,-split])
    
    # Estimate est_x_miss for df_k_fold_miss
    mu_miss <- mu[split]
    mu_obs <- mu[-split]
    
    sigma_miss_obs <- sigma_trunc[split, -split]
    sigma_obs_obs <- sigma_trunc[-split, -split]
    
    est_x_miss <- lapply(1:l1, function(n){mu_miss+sigma_miss_obs%*%solve(sigma_obs_obs)%*%(df_k_fold_obs[n,]-mu_obs)})
    
    # Error of estimated and true missing observation
    for(s in 1:l1){
      mse1[k] <- mse1[k] + norm(est_x_miss[[s]]-df_k_fold_miss[s,], type = "2")^2/l1/K
    }
  }
  
  return(sum(mse1))
}

WrongPCAImproved(df1, 1:2, 5)


# Cross-Validation: "Missing Data Approach"

# Split data into missing and observed data on jth component
# e.g. split missing and observed data in half, or split variables randomly
# splitting the data can be done for all variables regardless of r
# split <- sample(1:p, ceiling(ratio)), ratio !< r
split <- 1:ceiling(r/2)

mse2 <- rep(0, K)

for (k in 1:K) {
  l1 <- length(samples[[k]])
  
  df_k <- df[-samples[[k]],]
  mu <- colMeans(df_k)
  svd_sigma <- svd(cov(df_k))
  sigma_trunc <- svd_sigma$u[,-trunc] %*% diag(svd_sigma$d[-trunc]) %*% t(svd_sigma$u[,-trunc])

  df_k_fold <- df[samples[[k]],]
  df_k_fold_miss <- as.matrix(df_k_fold[,split])
  df_k_fold_obs <- as.matrix(df_k_fold[,-split])
  
  # Estimate est_x_miss for df_k_fold_miss
  mu_miss <- mu[split]
  mu_obs <- mu[-split]
  sigma_miss_obs <- sigma_trunc[split, -split]
  sigma_obs_obs <- sigma_trunc[-split, -split]
  
  est_x_miss <- lapply(1:l1, function(n){mu_miss+sigma_miss_obs%*%solve(sigma_obs_obs)%*%(df_k_fold_obs[n,]-mu_obs)})
  
  # Error of estimated and true missing observation
  for(s in 1:l1){
    mse2[k] <- mse2[k] + norm(est_x_miss[[s]]-df_k_fold_miss[s,], type = "2")^2/l1/K
  }
  
}


KDEApproach <- function(X, pca_vec){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r)
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  sigma <- cov(X)
  svd_sigma <- svd(sigma)
  sigma_trunc <- svd_sigma$u[,pca_vec] %*% diag(svd_sigma$d[pca_vec]) %*% t(svd_sigma$v)[pca_vec,]
  
  return(norm(sigma-sigma_trunc, type = "F"))
}

KDEApproach(df1, 1:2)


MatrixCompletion <- function(X, pca_vec){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r)
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  # Create bivariate index set for observed data
  n <- dim(X)[1]
  p <- dim(X)[2]
  biv_index <- matrix(sample(0:1, size = n*p, replace = T), nrow = n, ncol = p)
  
  # Iterative-hard thresholding algorithm
  M <- X
  for (i in 1:10) {
    svd_M <- svd(M)
    M_trunc <- svd_M$u[,pca_vec] %*% diag(svd_M$d[pca_vec]) %*% t(svd_M$v)[pca_vec,]
    M <- X * biv_index + M_trunc * !biv_index
  }
  
  return(norm(X-M, type = "F"))
}
  
MatrixCompletion(df1, 1:2)




# Question:

# 25.11.22
# How is r choosen? -> PCA first
# How missing and observed data split? -> randomly, but rule of thumb is half half
# MSE on x_miss? -> Yes
# Less informations in Method 1? -> No, different computations of mu and sigma
# Method 1: Truncating Missed and Observed data accordingly? -> No
# Method 2: Still p dimensions? -> Yes
# Method 3: How to cross-validate? Different error measure of covariance similiarity

# ..
# Method 0: Correctly implemented?
# Method 1: Splitting missing and observed data differently in Folds?
# Method 4: When is convergence achieved? How to choose tolerance?










