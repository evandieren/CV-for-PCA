# mvtnorm provides possibility to create multivariate gaussian rv
library(mvtnorm)
library(MASS)
library(psych)
library(matrixStats)
library(expss)

EM_MissingData <- function(X,matrix_miss){
  n <- nrow(X)
  p <- ncol(X)
  data <- X
  
  for (i in 1:n){
    data[i,matrix_miss[i,]] <- NA
  }
  
  # Initialize the parameters (Initial guesses)
  mu_hat <- colMeans(data, na.rm = TRUE)
  sigma_hat <-  diag(p)
  
  X_hat <- data
  sigma_hat_old <- sigma_hat+600*diag(p)
  tol <- 0.01
  while (norm(sigma_hat-sigma_hat_old,type="F")/norm(sigma_hat_old,type="F") > tol){
    #print(norm(sigma_hat-sigma_hat_old,type="F")/norm(sigma_hat_old,type="F"))
    sigma_hat_old <- sigma_hat
    
    # Filling the missing values with the linear methods (m_i)
    for (i in 1:n){
      index_miss <- matrix_miss[i,]
      X_hat[i,index_miss] <- mu_hat[index_miss] + sigma_hat[index_miss,-index_miss]%*%ginv(sigma_hat[-index_miss,-index_miss])%*%(X_hat[i,-index_miss]-mu_hat[-index_miss])
    }
    
    mu_hat <- colMeans(X_hat)
    sigma_hat <- cov(X_hat)
  }
  out = list("mu"=mu_hat,"sig"=sigma_hat)
  return(out)
}


EM_MissingData_Monica <- function(X,matrix_miss,n_iterations){
  n <- nrow(X)
  p <- ncol(X)

  data <- X
  
  for (i in 1:n){
    data[i,matrix_miss[i,]] <- NA
  }
  
  
  # Initialize the parameters
  mu_hat <- colMeans(data, na.rm = TRUE)
  sigma_hat <-  diag(p) #cov(data, use = "complete.obs")
  
  #mu_hat <- rep(0, p)
  #mat <- matrix(rnorm(p*p, mean = 0, sd = ), nrow = p, ncol = p)
  #sigma_hat <- diag(p) #mat %*% t(mat) 
  X_hat <- data
  sigma_hat_old <- sigma_hat+600*diag(p)
  tol <- 0.01
  while (norm(sigma_hat-sigma_hat_old,type="F")/norm(sigma_hat_old,type="F") > tol){
    print(norm(sigma_hat-sigma_hat_old,type="F")/norm(sigma_hat_old,type="F"))
    
    sigma_hat_old <- sigma_hat
    # Filling the missing values with the linear methods (m_i)
    #first_term <- matrix(0,p,p)
    for (i in 1:n){
      index_miss <- matrix_miss[i,]
      X_hat[i,index_miss] <- mu_hat[index_miss] + sigma_hat[index_miss,-index_miss]%*%ginv(sigma_hat[-index_miss,-index_miss])%*%(X_hat[i,-index_miss]-mu_hat[-index_miss])
      #V_i <- sigma_hat[index_miss,index_miss]-sigma_hat[index_miss,-index_miss]%*%ginv(sigma_hat[-index_miss,-index_miss])%*%sigma_hat[-index_miss,index_miss]
      #C_top_l <- X_hat[i,index_miss]%*%t(X_hat[i,index_miss])+V_i
      #print(dim(C_top_l))
      #C_top_r <- X_hat[i,index_miss]%*%t(X_hat[i,-index_miss])
      #print(dim(C_top_r))
      #C_low_l <- t(C_top_r)
      #print(dim(C_low_l))
      #C_low_r <- X_hat[i,-index_miss]%*%t(X_hat[i,-index_miss])
      #print(dim(C_low_r))
      #first_term <- first_term + rbind(cbind(C_top_l,C_top_r),cbind(C_low_l,C_low_r))/n
    } # X_obs et X_hat_miss ok -> E[x_i] est ok : it is x_hat
      mu_hat <- colMeans(X_hat)
      sigma_hat <- matrix(0,p,p)
      for (i in 1:n){
        sigma_hat <- sigma_hat + (X_hat[i,]-mu_hat)%*%t(X_hat[i,]-mu_hat)
      }
      #sigma_hat <- sigma_hat/(n-1)
      sigma_hat <- first_term - mu_hat%*%t(mu_hat)
    }
  out = list("mu"=mu_hat,"sig"=sigma_hat)
  return(out)
}





MissingData <- function(X, folds){
  # args: X matrix containing data, folds is the indexes of the rows contained in each fold
  # returns: MSE of the CV

  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(folds)

  matrix_miss <- matrix(data=0,n,floor(p/2))
  for (j in 1:n){
    matrix_miss[j,] <- sample(1:p, floor(p/2))
  }
    
  mse2 <- rep(0, p)
  for (r in 1:p) {
    for (k in 1:K) {
      
      n_fold <- length(folds[,k])
      df_k <- X[-folds[,k],]
      
      ls <- EM_MissingData(df_k, matrix_miss[-folds[,k],])
      mu <- ls$mu
      sigma <- ls$sig
      
      eigen_sigma <- eigen(sigma)
      sum_eigen <- sum(eigen_sigma$values[-c(1:r)]) # Eigenvalues destroyed
      eigen_sigma$values[-c(1:r)] <- 0 # truncation of the p-r last eigenvalues
      eigen_sigma_trunc <- eigen_sigma$vectors %*% diag(eigen_sigma$values) %*% t(eigen_sigma$vectors) + sum_eigen/p*diag(p)

      
      for (i in folds[,k]){ # looping on each row of the fold K
        index_miss <- matrix_miss[i,]
        index_obs <- (1:p)[-index_miss]
        mu_miss <- mu[index_miss]
        mu_obs <- mu[index_obs]
        sigma_miss_obs <- eigen_sigma_trunc[index_miss, index_obs]
        sigma_obs_obs <- eigen_sigma_trunc[index_obs, index_obs]
        est_x_n_miss <- mu_miss+sigma_miss_obs%*%ginv(sigma_obs_obs)%*%(X[i,index_obs]-mu_obs)
        
        mse2[r] <- mse2[r] + sum((est_x_n_miss-X[i,index_miss])^2)/n_fold
      }
    }
    mse2[r] <- mse2[r]/K
  }
  return(mse2)
}


# n <- 1000
# p <- 20
# r <- 7
# K <- 5
# X <- array(rnorm(n*p,0,1),c(n,p))
# svd_X <- svd(X)
# r_val <- svd_X$d[r]
# svd_X$d[-(1:r)] <- 0
# X_tronc <- svd_X$u %*% diag(svd_X$d) %*% t(svd_X$v) + rmvnorm(n = n, mean = rep(0,p), sigma = 0.001*r_val*diag(p))
# samples <- matrix(sample(1:n),ncol=K)
# 
# percent_miss = floor(p/2)
# matrix_miss <- matrix(data=0,n,percent_miss)
# for (j in 1:n)
#   matrix_miss[j,] <- sample(1:p,percent_miss)
# 
# out <- EM_MissingData(X_tronc,matrix_miss)
# print(norm(out$sig-cov(X_tronc),type="F")/norm(cov(X_tronc),type="F"))
# 
# out_miss <- MissingData(X_tronc, samples)
# 
# plot(1:p,out_miss,type="l")
# 
# print(which.min(out_miss))



EM_MissingDataProf <- function(X, matrix_miss){
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
  #l_comp <- Inf
  a <- +n/2*log(det(ginv))
  #l_comp_next <- a - 1/2*sum(sapply(1:n, function(m){tr(outer(X[m,]-mu, X[m,]-mu)%*%ginv)}))
  
  #while(abs(l_comp_next-l_comp)>tol){
  for (j in 1:50) {
    # Update observations and compute Pseudo-inverse of estimated sigma
    ginv <- ginv(sigma, tol=1e-20)
    for (i in 1:n){
      index_miss <- matrix_miss[i,] # B izarre qu'il y ait des doublons
      X[i,index_miss] <- mu[index_miss] + sigma[index_miss,-index_miss]%*%ginv[-index_miss,-index_miss]%*%(X[i,-index_miss]-mu[-index_miss])
    }
    
    # M-step
    C <- cov(X)
    mu <- colMeans(X)
    sigma <- (Reduce("+", lapply(1:n, function(m){outer(X[m,]-mu, X[m,]-mu)})) + C)/n
    
    # Likelihood for convergence
    #l_comp <- l_comp_next
    #l_comp_next <- +n/2*log(det(ginv(sigma, tol=1e-20))) - 1/2*sum(sapply(1:n, function(m){tr(outer(X[m,], X[m,])%*%ginv)}))
  }
  return(list(mu, sigma))
}

