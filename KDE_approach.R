library(mvtnorm)
KDEApproach <- function(X, pca_vec){
  # args: X matrix containing data, pca_vec containing index of pca (dimension r)
  # returns: error measure of truncated covariance to the real covariance (Frobenius norm)
  
  sigma <- cov(X)
  svd_sigma <- svd(sigma)
  sigma_trunc <- svd_sigma$u[,pca_vec] %*% diag(svd_sigma$d[pca_vec]) %*% t(svd_sigma$v)[pca_vec,]
  
  return(norm(sigma-sigma_trunc, type = "F"))
}


n <- 1000
mean_df <- rep(0, 4)


sim <- 1000
a <- matrix(rep(0, sim*4), nrow=sim, ncol = 4)
b <- matrix(rep(0, sim*4), nrow=sim, ncol = 4)
c <- matrix(rep(0, sim*4), nrow=sim, ncol = 4)
d <- matrix(rep(0, sim*4), nrow=sim, ncol = 4)
e <- matrix(rep(0, sim*4), nrow=sim, ncol = 4)
#similar to what Noah did

for (r in 2:5){
  for (i in 1:sim) {
    K = 10
    pca_vec <- 1:r
    mat <- matrix(rnorm(r*r, mean = 0, sd = 5), nrow = r, ncol = r)
    df <- rmvnorm(n = n, mean = mean_df[1:r], sigma = mat %*% t(mat))
    df_uni_high_noise <- df + rmvnorm(n = n, mean = mean_df[1:r], sigma = diag(rep(150, r)))
    df_uni_low_noise <- df + rmvnorm(n = n, mean = mean_df[1:r], sigma = diag(rep(1.5, r)))
    dec <- runif(n = r, min = 0, max = 1000)
    df_diff_noise <- df + rmvnorm(n = n, mean = mean_df[1:r], sigma = diag(dec))
    dec <- c(1:r)^2
    df_incr_noise <- df + rmvnorm(n=n, mean=mean_df[1:r], sigma=diag(dec))
    
    a[i,r-1] <- KDEApproach(df_uni_high_noise, pca_vec)
    b[i,r-1] <- KDEApproach(df_uni_low_noise, pca_vec)
    c[i,r-1] <- KDEApproach(df_diff_noise, pca_vec)
    d[i,r-1] <- KDEApproach(df_incr_noise, pca_vec)
    e[i,r-1] <- KDEApproach(df, pca_vec)
  }
}

