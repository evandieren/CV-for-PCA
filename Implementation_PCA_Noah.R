# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)


# Dataframe as input
df <- c()


# Create folds for Cross-Validation

df_df <- data.frame("obs1"=c(5,43,2,3,1,4,8,4,5,9), "obs2"=c(5,6,7,3,9,9,12,1,5,8), "obs3"=c(1,2,5,3,4,5,6,1,5,7), "obs4"=c(1,7,6,3,4,5,6,1,7,10))
df <- as.matrix(df_df)
n <- dim(df)[1]
p <- dim(df)[2]
index <- c(1:n)
K <- 3
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

# Index can now be used to slice df in every step of cross-validation

# How is r set? r < p
r <- 2

# Sample randomly p-r columns in order to keep r principal columns
trunc <- sample(1:p, p-r)



# Cross-Validation: "Artificially turn the unsupervised problem into a supervised one"

# Split data into missing and observed data on jth component
# e.g. split missing and observed data in half, or split variables randomly
# split <- sample(1:p, ceiling(ratio)), ratio < r as data is truncated to rank r
split <- sample(1:p, ceiling(p/2))

mse1 <- rep(0, K)

for (k in 1:K) {
  l1 <- length(samples[[k]])
  
  df_k <- df[-samples[[k]],]
  mu <- colMeans(df_k)
  svd_sigma <- svd(cov(df_k))
  sigma_trunc <- svd_sigma$u[,-trunc] %*% diag(svd_sigma$d[-trunc]) %*% t(svd_sigma$v)[-trunc,]
  
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
    mse1[k] <- mse1[k] + norm(est_x_miss[[s]]-df_k_fold_miss[s,], type = "2")^2/l1/K
  }
}



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


# KDE Approach

trunc <- sample(1:p, p-r)

KDEApproach <- function(X, pca_vec){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r)
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  sigma <- cov(X)
  svd_sigma <- svd(sigma)
  sigma_trunc <- svd_sigma$u[,pca_vec] %*% diag(svd_sigma$d[pca_vec]) %*% t(svd_sigma$v)[pca_vec,]
  
  return(norm(sigma-sigma_trunc, type = "F"))
}

KDEApproach(df, 1:2)


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
  
MatrixCompletion(df, 1:2)




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
# Method: Splitting missing and observed data differently in Folds?










