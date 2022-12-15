# This R file contains the simulation study conducted on different data sets
# for the Methods Wrong CA and WrongPCAImproved implemented in the file
# Implementation_PCA
library(mvtnorm)
source("Implementation_PCA.R")

# Goal:Create multiple very different data sets to perform simulations on them
# Simulation study on cov matrix(rnorm(p*p, mean = 0, sd = 12), nrow = p, ncol = p)
# !!! Takes a while to compute, reduce n and sim for faster overview

n <- 200
p <- 5
K = 5
mean_df <- rep(0, p)

sim <- 2
lsEigen <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsWrongPCA <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsWrongPCAImproved <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsKDEApproach <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsMatrixCompletion <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)

for (i in 1:sim) {
  samples <- matrix(sample(1:n),ncol=K)
  
  # mat <- matrix(runif(p*p, min = -20, 20), nrow = p, ncol = p)
  mat <- matrix(rnorm(p*p, mean = 0, sd = 1), nrow = p, ncol = p)
  df <- rmvnorm(n = n, mean = mean_df, sigma = mat %*% t(mat))
  
  # Gaussian data with uniform "high noise", "low noise", differing noise, increasing noise
  df_uni_high_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(rep(50, p)))
  df_uni_low_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(rep(2.2, p)))
  df_diff_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(runif(n = p, min = 0, max = 10)))
  df_incr_noise <- df + rmvnorm(n=n, mean=mean_df, sigma=diag(c(1:p)))
  data <- list("High noise"=df_uni_high_noise, "Low noise"=df_uni_low_noise, "Differing noise"=df_diff_noise, "Increasing noise"=df_incr_noise)
  
  for (da in 1:length(data)) {
    sin <- svd(data[[da]]) # svd because of asymmetric data set
    lsEigen[[da]][i,] <- sin$d
    lsWrongPCA[[da]][i,] <- WrongPCA(data[[da]], samples)
    lsWrongPCAImproved[[da]][i,] <- WrongPCAImproved(data[[da]], samples)
    lsKDEApproach[[da]][i,] <- KDEApproach(data[[da]])
    lsMatrixCompletion[[da]][i,] <- MatrixCompletion(data[[da]])
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

for (i in 1:4) {plot(1:p, colMeans(lsMatrixCompletion[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data)[i]), type = "b", pch = 19, lty = 1, col = 1)}
title("Error for Matrix Completion", line = - 1, outer = TRUE)


