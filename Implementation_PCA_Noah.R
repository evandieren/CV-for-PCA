# Dataframe as input
df <- c()


# Create folds for Cross-Validation

df_df <- data.frame("obs1"=c(5,43,2,3), "obs2"=c(5,6,7,3))
df <- as.matrix(df_df)
n <- dim(df)[1]
index <- c(1:n)
K <- 2
samples <- vector("list", K)

for (i in 1:K) {
  if (i < K) {
    fold <- sample(index, n/K)
    samples[[i]] <- fold
    index <- index[-fold]
  }
  else{
    samples[[i]] <- index
  }
}

# Index can now be used to slice df in every step of cross-validation step


