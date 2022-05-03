library(RegTree)
library(randomForest)
library(ranger)
library(dateutils)
library(devtools)
load_all()

# ------ Simulate Data --------------
b <- c(1,5,3,0,2,-1,2,-3,0,1,1,0,2,0,3,-1,-2,-1,2,4)
SimData <- function(n = 100){
  X <- matrix(rnorm(n*20),n,20)
  y <- X%*%b + rnorm(n)
  return(list(X = X,
              b = b,
              y = y))
}

n <- 200
# Run a single simulation
sim <- SimData(n)
X_train <- sim$X
y_train <- sim$y
sim <- SimData(n)
X_fit <- sim$X
y_true <- sim$y

# Estimate a single tree
to_keep <- seq(0,(n-1)) # do not use "bagging" for single tree, i.e. keep all obs
Tree <- Reg_Tree(y_train, X_train, to_keep)
# Get out of sample fitted values and feature contributions from this tree:
#  (in C++ iterating over columns is much faster than iterating over rows)
out <- FitMat(t(X_train), Tree)
ts.plot(cbind(y_train, out[[1]]), col = c("black", "deeppink"), lwd=2) # look at fit (in sample so cheating)
# Look at feature contributions
FC <- t(out[[2]])
# Feature contributions add up 
head(rowSums(FC) + mean(y_train))
head(c(out[[1]]))

# compare this with the linear parameter values from the DGP:
comp <- (t(b)%x%matrix(1,n,1))*X_train

# comparing the first period (this is in sample):

tst <- cbind(comp[1,], FC[1,])
colnames(tst) <- c("true contribution", "estimated contribution")
tst

