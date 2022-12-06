<<<<<<< HEAD
library(mvtnorm)
# create data vectors
sample_size <- 1000                                       
sample_meanvector <- c(10, -2, 5, 14, 20)                                   
sample_covariance_matrix <- matrix(c(7, 4, 3, 2, 1, 4, 5, 4, 3, 2,
                                     3, 4, 5, 4, 3, 2, 3, 4, 5, 4, 1, 
                                     2, 3, 4, 9), ncol = 5)
# create multivariate normal distribution
sample_distribution <- rmvnorm(n = sample_size,
                               mean = sample_meanvector, 
                               sigma = sample_covariance_matrix)
# print distribution
head(sample_distribution)
#sample_distribution
#It has no sense to plot because we have 5 dimensions

=======
library(MASS)
# create data vectors
sample_size <- 1000                                       
sample_meanvector <- c(10, -2, 5, 10, 20)                                   
sample_covariance_matrix <- matrix(c(7, 4, 3, 2, 1, 4, 5, 4, 3, 2,
                                     3, 4, 5, 4, 3, 2, 3, 4, 5, 4, 1, 
                                     2, 3, 4, 9), ncol = 5)
# create multivariate normal distribution
sample_distribution <- mvrnorm(n = sample_size,
                               mu = sample_meanvector, 
                               Sigma = sample_covariance_matrix)
# print distribution
head(sample_distribution)
#sample_distribution
>>>>>>> 320d26ddce8e291c3956038ebe2c1213caf83d97
