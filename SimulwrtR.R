library(mvtnorm)
library(latex2exp)
source("Methods/WrongPCA.R")
source("Methods/WrongPCAImproved.R")
source("Methods/Missing_data.R")
source("Methods/KDEApproach.R")
source("Methods/MatrixCompletion.R")

SimulationStudy <- function(method, str, n, p, K, r, sim, noise, meth=T, eigen=F){
  # method: CV method choosen
  # n: number of observations of data set
  # p: number of variables of data set
  # K: number of folds for cross-validation
  # r: rank of truncated data set
  # sim: amount of simulation runs
  # noise: vector consisting on high noise, low noise, differing noise and increasing noise
  # returns:
  #   plots and acceptance rate of dimensions for different data sets
  
  # Store scree plots as lists of every data set
  if (eigen){
    lsEigen0 <- replicate(3, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
    lsEigen1 <- replicate(3, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
    lsEigen2 <- replicate(3, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
  }
  
  if (meth){
    # Store CV Errors as lists of every data set
    lsmethod0 <- replicate(3, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
    lsmethod1 <- replicate(3, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
    lsmethod2 <- replicate(3, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
    
    # Store choosen dimensions of every data set
    choose0 <- replicate(3, rep(0, p), simplify=F)
    choose1 <- replicate(3, rep(0, p), simplify=F)
    choose2 <- replicate(3, rep(0, p), simplify=F)
  }

  dataset <- function(df,r,hn,ln,dn,incn){
    df_svd <- svd(df)
    df_svd$d[-(1:r)] <- 0
    last_sv <- df_svd$d[r]
    df <- df_svd$u %*% diag(df_svd$d) %*% t(df_svd$v)
    
    df_uni_high_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = hn*last_sv*diag(p))
    df_uni_low_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = ln*last_sv*diag(p))
    df_diff_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = diag(runif(n = p, min = 0, max = dn*last_sv)))
    #df_incr_noise <- df + rmvnorm(n=n, mean=rep(0, p), sigma=diag(last_sv*incn*c(1:p)/p))
    return(list("High noise"=df_uni_high_noise, "Low noise"=df_uni_low_noise, "Differing noise"=df_diff_noise))#, "Increasing noise"=df_incr_noise))
  }
  
  for (i in 1:sim) {
    samples <- matrix(sample(1:n),ncol=K)
    
    df <- rmvnorm(n = n, mean = rep(0, p), sigma = diag(p)) # basis real data set 0
    data0 <- dataset(df,r,noise[1],noise[2],noise[3])#,noise[4])
    
    mat <- matrix(rnorm(p*p, mean = 0, sd = 1), nrow = p, ncol = p)
    df <- rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat)) # basis real data set 1
    data1 <- dataset(df,r,noise[1],noise[2],noise[3])#,noise[4])
    
    mat <- matrix(1:p^2, nrow = p, ncol = p)/p
    df <- rmvnorm(n = n, mean = rep(0, p), sigma = mat %*% t(mat)) # basis real data set 2
    data2 <- dataset(df,2,noise[1],noise[2],noise[3])#,noise[4])
    
    for (da in 1:3) {
      if(eigen){
        lsEigen0[[da]][i,] <- svd(data0[[da]])$d
        lsEigen1[[da]][i,] <- svd(data1[[da]])$d
        lsEigen2[[da]][i,] <- svd(data2[[da]])$d
      }
      
      if(meth){
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
  }
  
  
  # Scree plots
  if(eigen){
    png("Figures_SimStudy/Scree_plot_0.png", width=900)#, res=res), pointsize=ps)
    par(mfrow=c(1,3))
    for (i in 1:3) {plot(1:p, colMeans(lsEigen0[[i]]), xlab="kth Singular value", ylab="Value", cex.lab = 1.5, main=TeX(paste0("Scree plot $D^",i,"_0$"), bold=T), cex.main = 2, type = "b", pch = 19, lty = 1, col = 1)}
    # title(TeX("Scree plot $D_0$", bold=T), line = - .9, outer = TRUE)
    dev.off()
    
    png("Figures_SimStudy/Scree_plot_1.png", width=900)#, res=res), pointsize=ps)
    par(mfrow=c(1,3))
    for (i in 1:3) {plot(1:p, colMeans(lsEigen1[[i]]), xlab="kth Singular value", ylab="Value", cex.lab = 1.5, main=TeX(paste0("Scree plot $D^",i,"_1$"), bold=T), cex.main = 2, type = "b", pch = 19, lty = 1, col = 1)}
    # title(TeX("Scree plot $D_1$", bold=T), line = - .9, outer = TRUE)
    dev.off()
    
    png("Figures_SimStudy/Scree_plot_2.png", width=900)#, res=res), pointsize=ps)
    par(mfrow=c(1,3))
    for (i in 1:3) {plot(1:p, colMeans(lsEigen2[[i]]), xlab="kth Singular value", ylab="Value", cex.lab = 1.5, main=TeX(paste0("Scree plot $D^",i,"_2$"), bold=T), cex.main = 2, type = "b", pch = 19, lty = 1, col = 1)}
    # title(TeX("Scree plot $D_2$", bold=T), line = - .9, outer = TRUE)
    dev.off()
    
  }

  if (meth) {
    # png(paste0("Figures_SimStudy/", str, "0.png"), width=900)#, res=res), pointsize=ps)
    # plot(1:p, colMeans(lsmethod0[[i]]), xlab="Rank r", ylab="Value", main=TeX(paste0("Error on $D^",i,"_0$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)

    # # Error of CV Methods
    # par(mfrow=c(1,3))
    # for (i in 1:3) {
    #   #png(paste0("Figures_SimStudy/", str, "0-",i,".png"), width=300)#, res=res), pointsize=ps)
    #   plot(1:p, colMeans(lsmethod0[[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(str," on $D^",i,"_0$"), bold=T), type = "b")
    #   #dev.off()
    #   }
    # #title(main=paste0("Error of ",str), line = - .9, outer = TRUE,cex.main=3)
    # #dev.off()
    # 
    # #png(paste0("Figures_SimStudy/", str, "1.png"), width=900)#, res=res), pointsize=ps)
    # par(mfrow=c(1,3))
    # for (i in 1:3) {plot(1:p, colMeans(lsmethod1[[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(str," on $D^",i,"_1$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}
    # title(paste0("Error of ",str), line = - .9, outer = TRUE)
    # #dev.off()
    # 
    # #png(paste0("Figures_SimStudy/", str, "2.png"), width=900)#, res=res), pointsize=ps)
    # par(mfrow=c(1,3))
    # for (i in 1:3) {plot(1:p, colMeans(lsmethod2[[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(str," on $D^",i,"_2$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}
    # title(paste0("Error of ",str), line = - .9, outer = TRUE)
    # #dev.off()
    
    c <- c()
    d <- c()
    chosen0 <- data.frame(matrix(ncol=p, nrow=3))
    colnames(chosen0) <- sapply(1:p, function(i){append(c, paste0("Var",i))})
    rownames(chosen0) <- c("0: High noise", "0: Low noise", "0: Differing noise")#, "0: Increasing noise")
    chosen0[1:3, 1:p] <- t(sapply(1:3, function(i){cbind(d,choose0[[i]]/sim)}))
    
    c <- c()
    d <- c()
    chosen1 <- data.frame(matrix(ncol=p, nrow=3))
    colnames(chosen1) <- sapply(1:p, function(i){append(c, paste0("Var",i))})
    rownames(chosen1) <- c("1: High noise", "1: Low noise", "1: Differing noise")#, "1: Increasing noise")
    chosen1[1:3, 1:p] <- t(sapply(1:3, function(i){cbind(d,choose1[[i]]/sim)}))
    
    c <- c()
    d <- c()
    chosen2 <- data.frame(matrix(ncol=p, nrow=3))
    colnames(chosen2) <- sapply(1:p, function(i){append(c, paste0("Var",i))})
    rownames(chosen2) <- c("2: High noise", "2: Low noise", "2: Differing noise")#, "2: Increasing noise")
    chosen2[1:3, 1:p] <- t(sapply(1:3, function(i){cbind(d,choose2[[i]]/sim)}))
    
    return(list(chosen0, chosen1, chosen2, lsmethod0, lsmethod1, lsmethod2))
    
  }
}



n <- 100
p <- 8
K <- 5
r <- 3
sim <- 5
noise <- c(0.02,0.001,0.005)
name <- "Matrix Completion"

# Uncomment to simulate :-)
# WrongPCA, WrongPCAImproved, MissingData, MatrixCompletion, KDEApproach
#for (r in 1:p){
#  set.seed(1312)
#  chosen <- SimulationStudy(MatrixCompletion,name, n, p, K, r, sim, noise, meth=T, eigen=F)
#  save(chosen,file=paste0("./datasets_plots/r_analysis/",name,"_",r,".Rdata"))
#}

set.seed(1312)
# chosen <- SimulationStudy(MatrixCompletion,name, n, p, K, r, sim, noise, meth=T, eigen=F)


ls <- list("Wrong PCA Improved"=WrongPCAImproved)
chosen <- SimulationStudy(ls$`"Wrong PCA Improved"`, names(list), n, p, K, r, sim, noise, meth=T, eigen=F)
par(mfrow=c(3,3))
for (i in 1:3) {plot(1:p, colMeans(chosen[[4]][[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(name," on $D^",i,"_0$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}
for (i in 1:3) {plot(1:p, colMeans(chosen[[5]][[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(name," on $D^",i,"_1$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}
for (i in 1:3) {plot(1:p, colMeans(chosen[[6]][[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(name," on $D^",i,"_2$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}


n <- 100
p <- 8
K <- 5
sim <- 5
noise <- c(0.02,0.001,0.005)

df <- rmvnorm(n = n, mean = rep(0, p), sigma = diag(p))

svd_df <- svd(df)

last <- svd_df$d[8]

df_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = 0.001*last*diag(p))

qr(df_noise)$rank

name <- "Matrix Completion"
# n <- 100
# p <- 8
# K <- 5
# r <- 3
# sim <- 5
# noise <- c(0.02,0.001,0.005)
# name <- "Matrix Completion"
# 
# # Uncomment to simulate :-)
# # WrongPCA, WrongPCAImproved, MissingData, MatrixCompletion, KDEApproach
# #for (r in 1:p){
# #  set.seed(1312)
# #  chosen <- SimulationStudy(MatrixCompletion,name, n, p, K, r, sim, noise, meth=T, eigen=F)
# #  save(chosen,file=paste0("./datasets_plots/r_analysis/",name,"_",r,".Rdata"))
# #}
# 
# set.seed(1312)
# # chosen <- SimulationStudy(MatrixCompletion,name, n, p, K, r, sim, noise, meth=T, eigen=F)
# 
# name <- "Wrong PCA Improved"
# chosen <- SimulationStudy(WrongPCAImproved, name, n, p, K, r, sim, noise, meth=T, eigen=F)
# par(mfrow=c(3,3))
# for (i in 1:3) {plot(1:p, colMeans(chosen[[4]][[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(name," on $D^",i,"_0$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}
# for (i in 1:3) {plot(1:p, colMeans(chosen[[5]][[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(name," on $D^",i,"_1$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}
# for (i in 1:3) {plot(1:p, colMeans(chosen[[6]][[i]]), xlab="Rank r", ylab="Error", main=TeX(paste0(name," on $D^",i,"_2$"), bold=T), type = "b", pch = 19, lty = 1, col = 1)}

