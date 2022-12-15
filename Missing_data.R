# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)
library(psych)

EM_MissingData <- function(X, matrix_miss){
  # X is data without kth fold
  # matrix_miss rowwise index of values missing in every observations
  
  # Initial stuff
  X_old <- X
  n <- nrow(X)
  p <- ncol(X)
  mu <- rep(0, p)
  mat <- matrix(rnorm(p*p, mean = 0, sd = 1), nrow = p, ncol = p)
  sigma <- mat %*% t(mat)
  ginv <- ginv(sigma, tol=1e-20)
  
  #Calculate missing values
  for (i in 1:n){
    index_miss <- matrix_miss[i,]
    X[i,index_miss] <- mu[index_miss] + sigma[index_miss,-index_miss]%*%ginv[-index_miss,-index_miss]%*%(X[i,-index_miss]-mu[-index_miss])
  }
  
  # Convergence stuff
  tol <- 1e-3
  l_comp <- Inf
  a <- +n/2*log(det(ginv))
  l_comp_next <- a - 1/2*sum(sapply(1:n, function(m){tr(outer(X[m,]-mu, X[m,]-mu)%*%ginv)}))
  
  while(abs(l_comp_next-l_comp)>tol){
    # Update observations and compute Pseudo-inverse of estimated sigma
    ginv <- ginv(sigma, tol=1e-20)
    for (i in 1:n){
      index_miss <- matrix_miss[i,]
      X[i,index_miss] <- mu[index_miss] + sigma[index_miss,-index_miss]%*%ginv[-index_miss,-index_miss]%*%(X[i,-index_miss]-mu[-index_miss])
    }
    
    # M-step
    C <- cov(X)
    mu <- colMeans(X)
    sigma <- (Reduce("+", lapply(1:n, function(m){outer(X[m,]-mu, X[m,]-mu)})) + C)/n
    
    # Likelihood for convergence
    l_comp <- l_comp_next
    l_comp_next <- +n/2*log(det(ginv(sigma, tol=1e-20))) - 1/2*sum(sapply(1:n, function(m){tr(outer(X[m,], X[m,])%*%ginv)}))
  }
  return(list(mu, sigma))
}



MissingData <- function(X, folds){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), K amount of CV-Folds
  # returns: MSE of the CV

  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(folds)
  matrix_miss <- matrix(sample(1:p, n*floor(p/2), replace=T), nrow=n)
  
  mse2 <- rep(0, p)
  for (r in 1:p) {
    for (k in 1:K) {
      
      n_fold <- length(folds[,k])
      df_k <- X[-folds[,k],]
      
      ls <- EM_MissingData(df_k, matrix_miss[-folds[,k],])
      mu <- ls[[1]]
      sigma <- ls[[2]]
      
      eigen_sigma <- eigen(sigma)
      eigen_sigma$values[-c(1:r)] <- 0 # truncation of the p-r last eigenvalues
      eigen_sigma_trunc <- eigen_sigma$vectors %*% diag(eigen_sigma$values) %*% t(eigen_sigma$vectors)
      
      # split observations
      index_miss <- sample(1:p, floor(p/2))
      index_obs <- (1:p)[-index_miss]
      mu_miss <- mu[index_miss]
      mu_obs <- mu[index_obs]
      sigma_miss_obs <- eigen_sigma_trunc[index_miss, index_obs]
      sigma_obs_obs <- eigen_sigma_trunc[index_obs, index_obs]
      
      # predict missing values
      df_fold <- X[folds[,k],]
      est_x_miss <- lapply(1:n_fold, function(n){mu_miss+sigma_miss_obs%*%ginv(sigma_obs_obs, tol = 1e-20)%*%(df_fold[n,index_obs]-mu_obs)})
      mse2[r] <- sum(sapply(1:n_fold, function(s){norm(est_x_miss[[s]]-df_fold[s,index_miss], type = "2")^2/n_fold})) + mse2[r]
    }
  }
  return(mse2)
}

n <- 100
p <- 5
K <- 5
mat <- matrix(rnorm(p*p, mean = 0, sd = 12), nrow = p, ncol = p)
df1 <- as.matrix(rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat)))
folds <- matrix(sample(1:n),ncol=K)
matrix_miss <- matrix(sample(1:p, n*floor(p/2), replace=T), nrow=n)
debug(EM_MissingData)
a <- EM_MissingData(df1, matrix_miss)


MissingData(df1, folds)
