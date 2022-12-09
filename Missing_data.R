# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)
library(psych)

df_df <- data.frame("obs1"=c(5,43,2,3,1,4,8,4,5,9), "obs2"=c(5,43,2,3,1,4,8,4,5,9), "obs3"=c(1,2,5,3,4,5,6,1,5,7), "obs4"=c(1,7,6,3,4,5,6,1,7,10))
df <- as.matrix(df_df)


# Create Multivariate Random Variables

EM_MissingData <- function(X, matrix_miss){

  n <- nrow(X)
  p <- ncol(X)
  
  # Set inital
  mu <- rep(0, p)
  sigma <- matrix(rep(1, p*p), nrow=p)
  
  #Calculate missing values
  for (i in 1:n){
    index_miss <- matrix_miss[i,]
    X[i,index_miss] <- mu[index_miss] + sigma[index_miss,-index_miss]%*%ginv(-index_miss,-index_miss, tol=1e-20)%*%(mu[-index_miss]-X[-index_miss])
  }
  
  tol <- 1e-5
  l_comp <- Inf
  # maybe modify to use tr()
  l_comp_next <- -N/2*log(det(ginv(sigma, tol=1e-20))) - 1/2*sum(sapply(1:n, function(m){tr(outer(X[m,], X[m,])%*%ginv(sigma, tol=1e-20))}))
  while(abs(l_comp_next - l_comp)>tol){
    # Utility variable
    var <- sum(sapply(1:n, function(m){tr(outer(X[m,], X[m,])%*%ginv(sigma, tol=1e-20))}))
    outer <- lapply(1:n, )
    
    # estimate missing values
    for (i in 1:n){
      index_miss <- matrix_miss[i,]
      X[i,index_miss] <- mu[index_miss] + sigma[index_miss,-index_miss]%*%ginv(-index_miss,-index_miss, tol=1e-20)%*%(mu[-index_miss]-X[-index_miss])
    }
    # Compute E-Step
    Q <- -N/2*log(det(ginv(sigma, tol=1e-20))) - 1/2*var
    
    # Update parameters
    mu <- colMeans(X)
    sigma <- Reduce("+", lapply(1:n, function(m){outer(X[m,]-mu, X[m,]-mu)})) + 
    
    # Calculate likelihood for convergence
    l_comp <- l_comp_next
    l_comp_next <- -N/2*log(det(ginv(sigma, tol=1e-20))) - 1/2*sum(sapply(1:n, function(m){tr(outer(X[m,], X[m,])%*%ginv(sigma, tol=1e-20))}))
    
  }
  return(list(mu, sigma))
}


MissingData <- function(X, folds){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), K amount of CV-Folds
  # returns: MSE of the CV

  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(folds)
  #index_miss <- sample(1:p, floor(p/2))
  #index_obs <- (1:p)[-index_miss]
  
  matrix_miss < matrix(sample(1:p, n*floor(p/2), replace=T), nrow=n)
  
  mse2 <- rep(0, p)
  for (r in 1:p) {
    for (k in 1:K) {

      n_fold <- length(folds[,k])
      df_k <- X[-folds[,k],]
      
      # EM-Algorithm to find sigma
      
      
      mu <- colMeans(df_k)
      eigen_sigma <- eigen(cov(df_k))
      eigen_sigma$values[-c(1:r)] <- 0 # truncation of the p-r last eigenvalues
      eigen_sigma_trunc <- eigen_sigma$vectors %*% diag(eigen_sigma$values) %*% t(eigen_sigma$vectors)
      
      df_fold <- X[folds[,k],]
      
      mu_miss <- mu[index_miss]
      mu_obs <- mu[index_obs]
      sigma_miss_obs <- eigen_sigma_trunc[index_miss, index_obs]
      sigma_obs_obs <- eigen_sigma_trunc[index_obs, index_obs]
      
      # Set tol = 1e-20 to reduce computational issues manually, don't know though if we allow sth that shouldn't be permitted
      est_x_miss <- lapply(1:n_fold, function(n){mu_miss+sigma_miss_obs%*%ginv(sigma_obs_obs, tol = 1e-20)%*%(df_fold[n,index_obs]-mu_obs)})
      
      # est_x_miss is a list so substracting from a matrix won't be possible
      # mse2[k] <- sum((est_x_miss-df_fold[,index_miss])^2)/(n_fold*dim_x) # Normalize with dimension?
      mse2[r] <- sum(sapply(1:n_fold, function(s){norm(est_x_miss[[s]]-df_fold[s,index_miss], type = "2")^2/n_fold})) + mse2[r]
    }
  }
  
  return(mse2)
}

n <- 1000
p <- 15
K <- 10
mat <- matrix(rnorm(p*p, mean = 0, sd = 12), nrow = p, ncol = p)
df1 <- rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat))
df1 <- as.matrix(df1)
folds <- matrix(sample(1:n),ncol=K)
MissingData(df1, folds)
a <- list(mat, df1)
