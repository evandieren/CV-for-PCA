library(mvtnorm)
source("WrongPCA.R")
source("WrongPCAImproved.R")
source("Missing_data.R")
source("KDEApproach.R")
source("MatrixCompletion.R")

dataset <- function(df,r,noi){
  df_svd <- svd(df)
  df_svd$d[-(1:r)] <- 0
  last_sv <- df_svd$d[r]
  df <- df_svd$u %*% diag(df_svd$d) %*% t(df_svd$v)
  
  df_uni_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = noi*last_sv*diag(p))
  return(list("Noise"=df_uni_noise))
}

SimulationStudy <- function(method, n, p, K, r, sim, noi){
  # method: CV method choosen
  # r: rank of truncated data set
  # sim: amount of simulation runs
  # return: plots and acceptance rate of dimensions for different data sets
  
  # Store scree plots as lists of every data set
  lsEigen0 <- replicate(1, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
 
  # Store CV Errors as lists of every data set
  lsmethod0 <- replicate(1, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  
  # Store choosen dimensions of every data set
  choose0 <- replicate(1, rep(0, p), simplify=F)
  
  for (i in 1:sim) {
    samples <- matrix(sample(1:n),ncol=K)
    
    df <- rmvnorm(n = n, mean = rep(0, p), sigma = diag(p)) # basis real data set 0
    data0 <- dataset(df,r,noi)
    
    for (da in 1:1) {
      lsEigen0[[da]][i,] <- svd(data0[[da]])$d
      
      lsmethod0[[da]][i,] <- method(data0[[da]], samples)
      
      min_err0 <- which.min(lsmethod0[[da]][i,])
      choose0[[da]][min_err0] <- choose0[[da]][min_err0] +1
      
    }
  }
  
  par(mfrow=c(1,2))
  # Scree plots
  plot1 = plot(1:p, colMeans(lsEigen0[[1]]), xlab="kth Eigenvalue", ylab="Value", main="Scree plot data set ", type = "b", pch = 19, lty = 1, col = 1)
  
  # Error of CV Methods
  plot2 = plot(1:p, colMeans(lsmethod0[[1]]), xlab="Rank r", ylab="Value", main= "Error data set ", type = "b", pch = 19, lty = 1, col = 1)
  
  c <- c()
  d <- c()
  chosen0 <- data.frame(matrix(ncol=p, nrow=1))
  colnames(chosen0) <- sapply(1:p, function(i){append(c, paste0("Var",i))})
  rownames(chosen0) <- c("Percentage")
  chosen0[1:1, 1:p] <- t(sapply(1:1, function(i){cbind(d,choose0[[i]]/sim)}))
 
  return(list(chosen0, lsmethod0))
}


n <- 100
p <- 7
K <- 5
sim <- 10
r <- 3 # truncated dimension of data set
noi <- 0.05

set.seed(1312)
# WrongPCA, WrongPCAImproved, MissingData, MatrixCompletion, KDEApproach
chosen <- SimulationStudy(WrongPCAImproved, n, p, K, r, sim, noi )

chosen[1]


