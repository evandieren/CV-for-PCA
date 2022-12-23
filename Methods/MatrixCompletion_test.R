source("MatrixCompletion.R")


data("Boston")

X <- as.matrix(Boston)[1:20, ]
nCV <- 14

set.seed(9302)
ind <- matrix(sample(1:length(X)), ncol = nCV)

# only first CV iteration
X_ <- X
X_[ind[,1]] <- NA

M <- lapply(1:ncol(X), function(r) MatrixCompletion_only(X = X_, X0 = X, r = r))

par(mfrow = c(2,2))
image(t(X_), main = "Incomplete matrix",)
image(t(M[[5]]), main = "Completed matrix")
image(t(X), main = "Original matrix")

# CV error on ommitted values
errors <- sapply(M, function(m) mean(sum((m[is.na(X_)]- X[is.na(X_)]))^2))
plot(errors, t = "b")


# Try with modified starting matrix ---------------------------------------
# Starting with M = X favours over-fitting, since the matrix completion starts
# at a unnaturally favorable initial variable.
# Start, e.g., at matrix imputing mean values for each variable instead:

X0 <- X
X0[ind[,1]] <- matrix(colMeans(X_, na.rm = T), byrow = TRUE, 
                      nrow = nrow(X), ncol = ncol(X))[ind[,1]]



M <- lapply(1:ncol(X), function(r) MatrixCompletion_only(X = X_, X0 = X0, r = r))

par(mfrow = c(2,2))
image(t(X_), main = "Incomplete matrix",)
image(t(M[[5]]), main = "Completed matrix")
image(t(X), main = "Original matrix")

# CV error on ommitted values
errors <- sapply(M, function(m) mean(sum((m[is.na(X_)]- X[is.na(X_)]))^2))
plot(errors, t = "b")
which.min(errors)
