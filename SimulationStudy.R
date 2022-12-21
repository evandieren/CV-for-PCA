library(mvtnorm)
source("WrongPCA.R")
source("WrongPCAImproved.R")
source("Missing_data.R")
source("KDEApproach.R")
source("MatrixCompletion.R")

SimulationStudy <- function(method, n, p, K, r, sim){
  # method: CV method choosen
  # r: rank of truncated data set
  # sim: amount of simulation runs
  # return: plots and acceptance rate of dimensions for different data sets
  
  # Store scree plots as lists of every data set
  lsEigen0 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  lsEigen1 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  lsEigen2 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  
  # Store CV Errors as lists of every data set
  lsmethod0 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  lsmethod1 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  lsmethod2 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  
  # Store choosen dimensions of every data set
  choose0 <- replicate(4, rep(0, p), simplify=F)
  choose1 <- replicate(4, rep(0, p), simplify=F)
  choose2 <- replicate(4, rep(0, p), simplify=F)
  
  dataset <- function(df,r,hn,ln,dn,incn){
    df_svd <- svd(df)
    df_svd$d[-(1:r)] <- 0
    last_sv <- df_svd$d[r]
    df <- df_svd$u %*% diag(df_svd$d) %*% t(df_svd$v)
    
    df_uni_high_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = hn*last_sv*diag(p))
    df_uni_low_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = ln*last_sv*diag(p))
    df_diff_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = diag(runif(n = p, min = 0, max = dn*last_sv)))
    df_incr_noise <- df + rmvnorm(n=n, mean=rep(0, p), sigma=diag(last_sv*incn*c(1:p)/p))
    return(list("High noise"=df_uni_high_noise, "Low noise"=df_uni_low_noise, "Differing noise"=df_diff_noise, "Increasing noise"=df_incr_noise))
  }
  
  for (i in 1:sim) {
    samples <- matrix(sample(1:n),ncol=K)
    
    df <- rmvnorm(n = n, mean = rep(0, p), sigma = diag(p)) # basis real data set 0
    data0 <- dataset(df,r,0,0,0,0)
    
    mat <- matrix(rnorm(p*p, mean = 0, sd = 1), nrow = p, ncol = p)
    df <- rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat)) # basis real data set 1
    data1 <- dataset(df,r,0.05,0.002,0.001,0.001)
    
    mat <- matrix(1:p^2, nrow = p, ncol = p)/p
    df <- rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat)) # basis real data set 2
    data2 <- dataset(df,2,0.05,0.002,0.001,0.001)
    
    for (da in 1:4) {
      lsEigen0[[da]][i,] <- svd(data0[[da]])$d
      lsEigen1[[da]][i,] <- svd(data1[[da]])$d
      lsEigen2[[da]][i,] <- svd(data2[[da]])$d
      
      lsmethod0[[da]][i,] <- method(data0[[da]], samples)
      lsmethod1[[da]][i,] <- method(data1[[da]], samples)
      lsmethod2[[da]][i,] <- method(data2[[da]], samples)
      
      min_err0 <- which.min(lsmethod0[[da]][i,])
      choose0[[da]][min_err0] <- choose0[[da]][min_err0] +1
      min_err1 <- which.min(lsmethod1[[da]][i,])
      choose1[[da]][min_err1] <- choose1[[da]][min_err1] +1
      min_err2 <- which.min(lsmethod2[[da]][i,])
      choose2[[da]][min_err2] <- choose2[[da]][min_err2] +1
    }
  }
  
  
  # Scree plots
  par(mfrow=c(1,4))
  for (i in 1:4) {plot(1:p, colMeans(lsEigen0[[i]]), xlab="kth Eigenvalue", ylab="Value", main=names(data0)[i], type = "b", pch = 19, lty = 1, col = 1)}
  title("Scree plot data set 0", line = - 1, outer = TRUE)
  for (i in 1:4) {plot(1:p, colMeans(lsEigen1[[i]]), xlab="kth Eigenvalue", ylab="Value", main=names(data1)[i], type = "b", pch = 19, lty = 1, col = 1)}
  title("Scree plot data set 1", line = - 1, outer = TRUE)
  for (i in 1:4) {plot(1:p, colMeans(lsEigen2[[i]]), xlab="kth Eigenvalue", ylab="Value", main=names(data2)[i], type = "b", pch = 19, lty = 1, col = 1)}
  title("Scree plot data set 2", line = - 1, outer = TRUE)
  
  # Error of CV Methods
  par(mfrow=c(1,4))
  for (i in 1:4) {plot(1:p, colMeans(lsmethod0[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data0)[i]), type = "b", pch = 19, lty = 1, col = 1)}
  title("Error data set 0", line = - 1, outer = TRUE)
  for (i in 1:4) {plot(1:p, colMeans(lsmethod1[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data1)[i]), type = "b", pch = 19, lty = 1, col = 1)}
  title("Error data set 1", line = - 1, outer = TRUE)
  for (i in 1:4) {plot(1:p, colMeans(lsmethod2[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data2)[i]), type = "b", pch = 19, lty = 1, col = 1)}
  title("Error data set 2", line = - 1, outer = TRUE)
  
  c <- c()
  d <- c()
  chosen0 <- data.frame(matrix(ncol=p, nrow=4))
  colnames(chosen0) <- sapply(1:p, function(i){append(c, paste0("Var",i))})
  rownames(chosen0) <- c("0: High noise", "0: Low noise", "0: Differing noise", "0: Increasing noise")
  chosen0[1:4, 1:p] <- t(sapply(1:4, function(i){cbind(d,choose0[[i]]/sim)}))
  
  c <- c()
  d <- c()
  chosen1 <- data.frame(matrix(ncol=p, nrow=4))
  colnames(chosen1) <- sapply(1:p, function(i){append(c, paste0("Var",i))})
  rownames(chosen1) <- c("1: High noise", "1: Low noise", "1: Differing noise", "1: Increasing noise")
  chosen1[1:4, 1:p] <- t(sapply(1:4, function(i){cbind(d,choose1[[i]]/sim)}))
  
  c <- c()
  d <- c()
  chosen2 <- data.frame(matrix(ncol=p, nrow=4))
  colnames(chosen2) <- sapply(1:p, function(i){append(c, paste0("Var",i))})
  rownames(chosen2) <- c("2: High noise", "2: Low noise", "2: Differing noise", "2: Increasing noise")
  chosen2[1:4, 1:p] <- t(sapply(1:4, function(i){cbind(d,choose2[[i]]/sim)}))
  
  return(list(chosen0, chosen1, chosen2, lsmethod0, lsmethod1, lsmethod2))
}


n <- 100
p <- 8
K <- 5
sim <- 10
r <- 3 # truncated dimension of data set

set.seed(1312)

# WrongPCA, WrongPCAImproved, MissingData, MatrixCompletion, KDEApproach
chosen <- SimulationStudy(KDEApproach, n, p, K, r, sim)
chosen[[1]]
chosen[[2]]
chosen[[3]]


