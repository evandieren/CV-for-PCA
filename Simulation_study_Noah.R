# This R file contains the simulation study conducted on different data sets
# for the Methods Wrong CA and WrongPCAImproved implemented in the file
# Implementation_PCA
library(mvtnorm)



set.seed(2)
# Goal:Create multiple very different data sets to perform simulations on them
n <- 1000
p <- 10
mean_df <- rep(0, p)

# mat <- matrix(runif(p*p, min = -20, 20), nrow = p, ncol = p)
mat <- matrix(rnorm(p*p, mean = 0, sd = 5), nrow = p, ncol = p)
# With mat computed as above: the following scree plots, i.e. look quite linear
# Hence pca-information of a balanced base matrix is hard to influence by different kind of noise
# Redundant information from scree plot and prcomp()

# mat <- matrix(1:p^2, nrow = p, ncol = p)/p^2
# Constant noise almost doesn't effect scree plot
# Differing noise strong uncharaterizable effects
# Increasing noise linearizes the scree plot
df <- rmvnorm(n = n, mean = mean_df, sigma = mat %*% t(mat))


# Gaussian data with uniform "high noise"
dec <- 20
uni_high_noise <- rmvnorm(n = n, mean = mean_df, sigma = diag(rep(dec, p)))
df_uni_high_noise <- df + uni_high_noise

# Scree plot
sin_high <- svd(df_uni_high_noise)
plot(1:p, sin_high$d)




# Gaussian data with uniform "low noise"
dec <- 2.2
uni_low_noise <- rmvnorm(n = n, mean = mean_df, sigma = diag(rep(dec, p)))
df_uni_low_noise <- df + uni_low_noise

# Scree plot
sin_low <- svd(df_uni_low_noise)
plot(1:p, sin_low$d)



# Gaussian data with differing noise in each component
dec <- runif(n = p, min = 0, max = 1000)
diff_noise <- rmvnorm(n = n, mean = mean_df, sigma = diag(dec))
df_diff_noise <- df + diff_noise

# Scree plot
sin_diff <- svd(df_diff_noise)
plot(1:p, sin_diff$d)



# Gaussian data with increasing noise
dec <- c(1:p)^2
diff_noise <- rmvnorm(n = n, mean = mean_df, sigma = diag(dec))
df_diff_noise <- df + diff_noise

# Scree plot
sin_diff <- svd(df_diff_noise)
plot(1:p, sin_diff$d)



# explicit test
p <- rmvnorm(40, mean = c(0, 0, 0), sigma = diag(c(1, 10, 13)))
plot(p[,1], p[,2], xlim = c(-15, 15), ylim = c(-15,15))
sin_p <- svd(p)
plot(1:3, sin_p$d)
pca1 <- prcomp(p)
summary(pca1)
p %*% pca1$rotation[,1:2]


# R built-in function for PCA
pca3 <- prcomp(df, rank. = p, tol = 0)
summary(pca3) # overview of the extent a variable explains the data
pca3$sdev # array sorted by Principal components that explain % of standard deviation




# Start simulation study, import functions from Implementation_PCA

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





