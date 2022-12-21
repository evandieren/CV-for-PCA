library(mvtnorm)
source("WrongPCAImproved.R")

n <- 1000
p <- 10
K = 10

sim <- 20

# Store scree plots as lists
lsEigen0 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsEigen1 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsEigen2 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)

# Store CV Errors as lists
lsWrongPCA0 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsWrongPCA1 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
lsWrongPCA2 <- replicate(4, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)

# Store choosen dimensions
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
  
  r <- 3
  df <- rmvnorm(n = n, mean = rep(0, p), sigma = diag(p)) # basis real data set 0
  # truncate to rank r

  data0 <- dataset(df,r,0.1,0.01,0.05,0.05)
  
  mat <- matrix(rnorm(p*p, mean = 0, sd = 1), nrow = p, ncol = p)
  df <- rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat)) # basis real data set 1
  data1 <- dataset(df,r,0.1,0.01,0.05,0.05)

  mat <- matrix(1:p^2, nrow = p, ncol = p)/p
  df <- rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat)) # basis real data set 2
  data2 <- dataset(df,2,0.1,0.01,0.05,0.05)
  
  for (da in 1:4) {
    lsEigen0[[da]][i,] <- svd(data0[[da]])$d
    lsEigen1[[da]][i,] <- svd(data1[[da]])$d
    lsEigen2[[da]][i,] <- svd(data2[[da]])$d
    
    lsWrongPCA0[[da]][i,] <- WrongPCA(data0[[da]], samples)
    lsWrongPCA1[[da]][i,] <- WrongPCA(data1[[da]], samples)
    lsWrongPCA2[[da]][i,] <- WrongPCA(data2[[da]], samples)
    
    min_err0 <- which(lsWrongPCA0[[da]][i,]==min(lsWrongPCA0[[da]][i,]))
    choose0[[da]][min_err0] <- choose0[[da]][min_err0] +1
    min_err1 <- which(lsWrongPCA1[[da]][i,]==min(lsWrongPCA1[[da]][i,]))
    choose1[[da]][min_err1] <- choose1[[da]][min_err1] +1
    min_err2 <- which(lsWrongPCA2[[da]][i,]==min(lsWrongPCA2[[da]][i,]))
    choose2[[da]][min_err2] <- choose2[[da]][min_err2] +1
  }
}


# Scree plots
par(mfrow=c(1,4))
for (i in 1:4) {plot(1:p, colMeans(lsEigen0[[i]]), xlab="kth Eigenvalue", ylab="Value", main=names(data)[i], type = "b", pch = 19, lty = 1, col = 1)}
title("Scree plot data set 0", line = - 1, outer = TRUE)
for (i in 1:4) {plot(1:p, colMeans(lsEigen1[[i]]), xlab="kth Eigenvalue", ylab="Value", main=names(data)[i], type = "b", pch = 19, lty = 1, col = 1)}
title("Scree plot data set 1", line = - 1, outer = TRUE)
for (i in 1:4) {plot(1:p, colMeans(lsEigen2[[i]]), xlab="kth Eigenvalue", ylab="Value", main=names(data)[i], type = "b", pch = 19, lty = 1, col = 1)}
title("Scree plot data set 2", line = - 1, outer = TRUE)

# Error of CV Methods
par(mfrow=c(1,4))
for (i in 1:4) {plot(1:p, colMeans(lsWrongPCA0[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data)[i]), type = "b", pch = 19, lty = 1, col = 1)}
title("Error data set 0", line = - 1, outer = TRUE)
for (i in 1:4) {plot(1:p, colMeans(lsWrongPCA1[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data)[i]), type = "b", pch = 19, lty = 1, col = 1)}
title("Error data set 1", line = - 1, outer = TRUE)
for (i in 1:4) {plot(1:p, colMeans(lsWrongPCA2[[i]]), xlab="Rank r", ylab="Value", main=paste0("Err ",names(data)[i]), type = "b", pch = 19, lty = 1, col = 1)}
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

test1 <- function(){
  print("test1")
}

test2 <- function(){
  print("test2")
}

usefunc <- function(x){
  x()
}

usefunc(test2)


