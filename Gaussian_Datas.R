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