# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)

df_df <- data.frame("obs1"=c(5,43,2,3,1,4,8,4,5,9), "obs2"=c(5,43,2,3,1,4,8,4,5,9), "obs3"=c(1,2,5,3,4,5,6,1,5,7), "obs4"=c(1,7,6,3,4,5,6,1,7,10))
df <- as.matrix(df_df)


# Create Multivariate Random Variables
sigma <- matrix(c(5,43,2,3,1,4,8,4,5,9,5,43,2,3,1,4,8,4,5,9,5,43,2,3,1), nrow = 5, ncol = 5)
df1 <- rmvnorm(n = 20, mean = rep(0, nrow(sigma)), sigma = sigma %*% t(sigma))

MissingData <- function(X, pca_vec, K){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r), K amount of CV-Folds
  # returns: MSE of the CV

  K <- K
  n_obs <- dim(X)[1]
  dim_x <- dim(X)[2]
  # pca_vec <- 1:r
  pca_vec <- 1:length(pca_vec)
  
  # Looks clean but I think observations can go missing, Maybe choosing ncol = K, is even smoother
  # maybe list where truncating sample(1:n_obs) on every between k, k+1 can also be implemented nicely
  folds <- matrix(sample(1:n_obs),ncol=n_obs/K)
  
  mse2 <- rep(0, K)
  for (k in 1:K) {
    # n_fold <- length(folds[k]) How does that work if folds is a Matrix? slice [k,] instead?
    n_fold <- length(folds[k,])
    df_k <- X[-folds[k,],]
    
    mu <- colMeans(df_k)
    eigen_sigma <- eigen(cov(df1))
    eigen_sigma$values[-pca_vec] = 0 # truncation of the p-r last eigenvalues
    eigen_sigma_trunc <- eigen_sigma$vectors %*% diag(eigen_sigma$values) %*% t(eigen_sigma$vectors)
    
    df_fold <- X[folds[k,],]
    index_miss <- sample(1:dim_x, floor(dim_x/2))
    index_obs <- (1:dim_x)[-index_miss]
    
    mu_miss <- mu[index_miss]
    mu_obs <- mu[index_obs]
    sigma_miss_obs <- eigen_sigma_trunc[index_miss, index_obs]
    sigma_obs_obs <- eigen_sigma_trunc[index_obs, index_obs]
    
    # Set tol = 1e-20 to reduce computational issues manually, don't know though if we allow sth that shouldn't be permitted
    est_x_miss <- lapply(1:n_fold, function(n){mu_miss+sigma_miss_obs%*%solve(sigma_obs_obs, tol = 1e-20)%*%(df_fold[n,index_obs]-mu_obs)})
    
    # est_x_miss is a list so substracting from a matrix won't be possible
    # mse2[k] <- sum((est_x_miss-df_fold[,index_miss])^2)/(n_fold*dim_x) # Normalize with dimension?
    mse2 <- sapply(1:n_fold, function(s){norm(est_x_miss[[s]]-df_fold[s,index_miss], type = "2")^2/(n_fold*dim_x)})
  }
  
  return(sum(mse2))
}

MissingData(df1, 1:2, 5)

