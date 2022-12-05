# This R file contains the simulation study conducted on different data sets
# for the Methods Wrong CA and WrongPCAImproved implemented in the file
# Implementation_PCA
library(mvtnorm)

# Goal:Create multiple very different data sets to perform simulations on them
n <- 1000
p <- 15
mean_df <- rep(0, p)

source("Implementation_PCA.R")

sim <- 1000
a <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
b <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
c <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
d <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
e <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
f <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
g <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
h <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
s1 <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
s2 <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
s3 <- matrix(rep(0, sim*p), nrow=sim, ncol = p)
s4 <- matrix(rep(0, sim*p), nrow=sim, ncol = p)

for (i in 1:sim) {
  K = 10
  samples <- matrix(sample(1:n),ncol=K)
  
  # mat <- matrix(runif(p*p, min = -20, 20), nrow = p, ncol = p)
  mat <- matrix(rnorm(p*p, mean = 0, sd = 12), nrow = p, ncol = p)
  # With mat computed as above: the following scree plots, i.e. look quite linearly decreasing
  # Hence pca-information of a balanced base matrix is hard to influence by different kind of noise
  # Redundant information from scree plot and prcomp()
  df <- rmvnorm(n = n, mean = mean_df, sigma = mat %*% t(mat))
  
  # Gaussian data with uniform "high noise", "low noise", differing noise, increasing noise
  df_uni_high_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(rep(50, p)))
  
  df_uni_low_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(rep(2.2, p)))
  
  dec <- runif(n = p, min = 0, max = 1000)
  df_diff_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(dec))
  
  dec <- c(1:p)^2
  df_incr_noise <- df + rmvnorm(n=n, mean=mean_df, sigma=diag(dec))
  
  sin_high <- svd(df_uni_high_noise)
  s1[i,] <- sin_high$d
  sin_low <- svd(df_uni_low_noise)
  s2[i,] <- sin_low$d
  sin_diff <- svd(df_diff_noise)
  s3[i,] <- sin_diff$d
  sin_incr <- svd(df_incr_noise)
  s4[i,] <- sin_incr$d
  
  a[i,] <- WrongPCA(df_uni_high_noise, samples)
  b[i,] <- WrongPCA(df_uni_low_noise, samples)
  c[i,] <- WrongPCA(df_diff_noise, samples)
  d[i,] <- WrongPCA(df_incr_noise, samples)
  
  e[i,] <- WrongPCAImproved(df_uni_high_noise, samples)
  f[i,] <- WrongPCAImproved(df_uni_low_noise, samples)
  g[i,] <- WrongPCAImproved(df_diff_noise, samples)
  h[i,] <- WrongPCAImproved(df_incr_noise, samples)
}



# Scree plots
par(mfrow=c(1,4))
plot(1:p, colMeans(s1), xlab="kth Eigenvalue", ylab="Value", main="High noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(s2), xlab="kth Eigenvalue", ylab="Value", main="Low noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(s3), xlab="kth Eigenvalue", ylab="Value", main="Differing noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(s4), xlab="kth Eigenvalue", ylab="Value", main="Increasing noise", type = "b", pch = 19, lty = 1, col = 1)
# mtext("Scree plots with different types of Gaussian noise", side=3, line=-1, outer=T, font="bold")
#dev.off()

par(mfrow=c(1,4))
plot(1:p, colMeans(a), xlab="Rank r", ylab="Value", main="Err high noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(b), xlab="Rank r", ylab="Value", main="Err low noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(c), xlab="Rank r", ylab="Value", main="Err differing noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(d), xlab="Rank r", ylab="Value", main="Err increasing noise", type = "b", pch = 19, lty = 1, col = 1)
# mtext("Scree plots with different types of Gaussian noise", side=3, line=-1, outer=T, font="bold")
#dev.off()

par(mfrow=c(1,4))
plot(1:p, colMeans(e), xlab="Rank r", ylab="Value", main="Err high noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(f), xlab="Rank r", ylab="Value", main="Err low noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(g), xlab="Rank r", ylab="Value", main="Err differing noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, colMeans(h), xlab="Rank r", ylab="Value", main="Err increasing noise", type = "b", pch = 19, lty = 1, col = 1)
# mtext("Scree plots with different types of Gaussian noise", side=3, line=-1, outer=T, font="bold")
#dev.off()




mat <- matrix(1:p^2, nrow = p, ncol = p)/p^2
# Constant noise almost doesn't effect scree plot
# Differing noise strong uncharaterizable effects
# Increasing noise linearizes the scree plot
df <- rmvnorm(n = n, mean = mean_df, sigma = mat %*% t(mat))

# Gaussian data with uniform "high noise", "low noise", differing noise, increasing noise
df_uni_high_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(rep(50, p)))

df_uni_low_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(rep(2.2, p)))

dec <- runif(n = p, min = 0, max = 1000)
df_diff_noise <- df + rmvnorm(n = n, mean = mean_df, sigma = diag(dec))

dec <- c(1:p)^2
df_incr_noise <- df + rmvnorm(n=n, mean=mean_df, sigma=diag(dec))

# Scree plots
par(mfrow=c(1,4))
sin_high <- svd(df_uni_high_noise)
plot(1:p, sin_high$d, xlab="kth Eigenvalue", ylab="Value", main="High noise", type = "b", pch = 19, lty = 1, col = 1)
sin_low <- svd(df_uni_low_noise)
plot(1:p, sin_low$d, xlab="kth Eigenvalue", ylab="Value", main="Low noise", type = "b", pch = 19, lty = 1, col = 1)
sin_diff <- svd(df_diff_noise)
plot(1:p, sin_diff$d, xlab="kth Eigenvalue", ylab="Value", main="Differing noise", type = "b", pch = 19, lty = 1, col = 1)
sin_incr <- svd(df_incr_noise)
plot(1:p, sin_incr$d, xlab="kth Eigenvalue", ylab="Value", main="Increasing noise", type = "b", pch = 19, lty = 1, col = 1)
# mtext("Scree plots with different types of Gaussian noise", side=3, line=-1, outer=T, font="bold")
#dev.off()

par(mfrow=c(1,4))
plot(1:p, WrongPCA(df_uni_high_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err high noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, WrongPCA(df_uni_low_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err low noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, WrongPCA(df_diff_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err differing noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, WrongPCA(df_incr_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err increasing noise", type = "b", pch = 19, lty = 1, col = 1)
# mtext("Scree plots with different types of Gaussian noise", side=3, line=-1, outer=T, font="bold")
#dev.off()

par(mfrow=c(1,4))
plot(1:p, WrongPCAImproved(df_uni_high_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err high noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, WrongPCAImproved(df_uni_low_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err low noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, WrongPCAImproved(df_diff_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err differing noise", type = "b", pch = 19, lty = 1, col = 1)
plot(1:p, WrongPCAImproved(df_incr_noise, samples), xlab="kth Eigenvalue", ylab="Value", main="Err increasing noise", type = "b", pch = 19, lty = 1, col = 1)
# mtext("Scree plots with different types of Gaussian noise", side=3, line=-1, outer=T, font="bold")
#dev.off()





