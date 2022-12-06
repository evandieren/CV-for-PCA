# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)


# df_df <- data.frame("obs1"=c(5,43,2,3,1,4,8,4,5,9), "obs2"=c(5,43,2,3,1,4,8,4,5,9), "obs3"=c(1,2,5,3,4,5,6,1,5,7), "obs4"=c(1,7,6,3,4,5,6,1,7,10))
# df <- as.matrix(df_df)


# Create Multivariate Random Variables
sigma <- matrix(c(5,43,2,3,1,4,8,4,5,9,5,43,2,3,1,4,8,4,5,9,5,43,2,3,1), nrow = 5, ncol = 5)
df1 <- rmvnorm(n = 20, mean = rep(0, nrow(sigma)), sigma = sigma %*% t(sigma))


# R built-in function for PCA
pca3 <- prcomp(df1, rank. = 5)
summary(pca3) # overview of the extent a variable explains the data
pca3$sdev # array sorted by Principal components that explain % of standard deviation


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
      a <- svd_X$u[,1:r]
      if (r==1){
        b <- as.matrix(svd_X$d[1:r])
      }
      else{
        b <- diag(svd_X$d[1:r])
      }
      c <- t(svd_X$v)[1:r,]
      X_trunc <- a %*% b %*% c
      U <- prcomp(X_trunc,rank.=r)$rotation
      
      # Define projection and project data
      P <- U%*%t(U)
      # X_trunc%*% solve(t(X_trunc) %*% X_trunc, tol = 1e-20) %*% t(X_trunc)
      df_k_proj <- P %*% t(df_k_fold)

      # Error of estimated and true missing observation
      mse[r] <- sum(sapply(1:l1, function(s){norm(df_k_fold[s,] - df_k_proj[,s], type = "2")^2/l1})) + mse[r]
      # Error is ridiculous
      # es <- df_k_proj
      # es1 <- df_k_fold
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
      l1 <- length(samples[,k])
      
      df_k <- X[-samples[,k],]
      mu <- colMeans(df_k)
      svd_sigma <- svd(cov(df_k))
      svd_sigma$d[-(1:r)] <- 0
      sigma_trunc <- svd_sigma$u %*% diag(svd_sigma$d) %*% t(svd_sigma$v)
      
      df_k_fold <- X[samples[,k],]
      df_k_fold_miss <- as.matrix(df_k_fold[,split])
      df_k_fold_obs <- as.matrix(df_k_fold[,-split])
      
      # Estimate est_x_miss for df_k_fold_miss
      mu_miss <- mu[split]
      mu_obs <- mu[-split]
      sigma_miss_obs <- sigma_trunc[split, -split]
      sigma_obs_obs <- sigma_trunc[-split, -split]
      
      est_x_miss <- lapply(1:l1, function(n){mu_miss + sigma_miss_obs %*% ginv(sigma_obs_obs, tol = 1e-20) %*% (df_k_fold_obs[n,]-mu_obs)})
      mse1[r] <- sum(sapply(1:l1, function(s){norm(est_x_miss[[s]]-df_k_fold_miss[s,], type = "2")^2/l1})) + mse1[r]
    }
  }
  return(mse1)
}


KDEApproach <- function(X, pca_vec){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r)
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  sigma <- cov(X)
  svd_sigma <- svd(sigma)
  sigma_trunc <- svd_sigma$u[,pca_vec] %*% diag(svd_sigma$d[pca_vec]) %*% t(svd_sigma$v)[pca_vec,]
  
  return(norm(sigma-sigma_trunc, type = "F"))
}


MatrixCompletion <- function(X, pca_vec){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r)
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  # Create bivariate index set for observed data
  n <- nrow(X)
  p <- ncol(X)
  biv_index <- matrix(sample(0:1, size = n*p, replace = T), nrow = n, ncol = p)
  
  # Iterative-hard thresholding algorithm
  M_old <- matrix(rep(1, size = n*p), nrow = n, ncol = p)
  M_l <- X
  tol = 1e-5
  while(norm(M_l-M_old, type = "F")>tol){
    #print("test")
    M_old = M_l
    svd_M <- svd(M_old)
    svd_M$d[-pca_vec] <- 0
    M_trunc <- svd_M$u %*% diag(svd_M$d) %*% t(svd_M$v)
    M_l <- X * biv_index + M_trunc * !biv_index
  }

  return(norm(X-M_l, type = "F"))
}



# Generate Folds for Cross-Validation

K <- 5
n <- nrow(df1)
index <- c(1:n)
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

#MatrixCompletion(df1, 1:1)

# samples <- matrix(sample(1:n_obs),ncol=n_obs/K)



# Test functions for Cross-Validation on same  Folds for consistency

# sin <- svd(df1) # Scree plot
# plot(1:5, sin$d)

# WrongPCA(df1, samples)
# WrongPCAImproved(df1, samples)

#WrongPCA_err <- sapply(2:5, function(k){WrongPCA(df1, samples)})
#WrongPCAImproved_err <- sapply(2:5, function(k){WrongPCAImproved(df1, samples)})
#plot(2:5, WrongPCA_err, "l", col = 1) # should return 0 for 5 but doesn't?
#lines(2:5, WrongPCAImproved_err, col = 2)

#KDEApproach_err <- sapply(2:5, function(n){KDEApproach(df1, 1:n)})
MatrixCompletion_err <- sapply(1:5, function(n){MatrixCompletion(df1, 1:n)})
#View(MatrixCompletion_err)
#plot(1:5, KDEApproach_err, "l", col = 1)
plot(1:5, log(MatrixCompletion_err),type="l", col = 2)

print(which.min(MatrixCompletion_err))

# Questions:

# 25.11.22
# How is r choosen? -> PCA first
# How missing and observed data split? -> randomly, but rule of thumb is half half
# MSE on x_miss? -> Yes
# Less informations in Method 1? -> No, different computations of mu and sigma
# Method 1: Truncating Missed and Observed data accordingly? -> No
# Method 2: Still p dimensions? -> Yes
# Method 3: How to cross-validate? Different error measure of covariance similiarity








