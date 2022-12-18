library(mvtnorm)
source("Implementation_PCA.R")
n <- 200
p <- 10
K = 5
mean_df <- rep(0, p)
sim = 2
lsEigen <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsWrongPCA <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsWrongPCAImproved <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsKDEApproach <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsMatrixCompletion <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)

for (i in 1:sim) {
  samples <- matrix(sample(1:n),ncol=K)
  mat <- diag(p)
  df <- rmvnorm(n = n, mean = mean_df, sigma = mat) #* t(matrix(rmvnorm(n = n, mean = mean_df, sigma = mat), nrow = p, ncol = n))
  
  df_0_25_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = (0.25*diag(p)))
  df_0_5_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = (0.5*diag(p)))
  df_0_75_noise <- df + rmvnorm(n = n, mean = mean_df, sigma =(0.75*diag(p)))
  df_1_noise <- df + rmvnorm(n = n, mean=mean_df, sigma= (diag(p)))
  data <- list("sigma 0.25 noise"=df_0_25_noise, "sigma 0.5 noise"=df_0_5_noise, "sigma 0.75 noise"=df_0_75_noise, "sigma 1 noise"=df_1_noise)

  for (da in 1:length(data)) {
    sin <- svd(data[[da]]) 
    lsEigen[[da]][i,] <- sin$d
    lsWrongPCA[[da]][i,] <- WrongPCA(data[[da]], samples)
    lsWrongPCAImproved[[da]][i,] <- WrongPCAImproved(data[[da]], samples)
    lsKDEApproach[[da]][i,] <- KDEApproach(data[[da]])
    #lsMatrixCompletion[[da]][i,] <- MatrixCompletion(data[[da]])
  }
}


# Scree plots
par(mfrow=c(1,4))
for (i in 1:4) {plot(1:p, colMeans(lsEigen[[i]]), xlab="kth Eigenvalue", ylab="Value", main=names(data)[i], type = "b", pch = 19, lty = 1, col = 1)}

# Error of CV Methods
par(mfrow=c(1,4))

for (i in 1:4) {plot(1:p, colMeans(lsWrongPCA[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data)[i]), type = "b", pch = 19, lty = 1, col = 1)}
title("Error for Wrong PCA", line = - 1, outer = TRUE)

for (i in 1:4) {plot(1:p, colMeans(lsWrongPCAImproved[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data)[i]), type = "b", pch = 19, lty = 1, col = 1)}
title("Error for Wrong PCA Improved", line = - 1, outer = TRUE)

for (i in 1:4) {plot(1:p, colMeans(lsKDEApproach[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data)[i]), type = "b", pch = 19, lty = 1, col = 1)}
title("Error for KDE Approach", line = - 1, outer = TRUE)
