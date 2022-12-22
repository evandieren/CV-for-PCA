percent_var <- function(Sigma){
  p <- dim(Sigma)[1]
  eigen_sig <- eigen(Sigma)
  V <- sum(eigen_sig$values)
  threshold <- 0.9
  for (r in 1:p){
    if (sum(eigen_sig$values[1:r]>threshold*V)){
      return(r)
    }
  }
}