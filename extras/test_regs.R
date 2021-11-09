x <- 1:40
y[1:20] <- x[1:20] + 2*rnorm(20) + 5
y[21:40] <- -x[21:40] + 4*rnorm(20) + 4

plot(x,y)

out <- find_cut(x, y, 1:40)
out[[1]]


bob <- quickreg(x[1:20],y[1:20],0)
bob[[1]]

# ------ Simulate Data --------------
SimData <- function(n = 100){
  X <- matrix(rnorm(n*20),n,20)
  b <- c(1,5,3,0,2,-1,2,-3,0,1,1,0,2,0,3,-1,-2,-1,2,4)
  y <- X%*%b + rnorm(n)
  return(list(X = X,
              b = b,
              y = y))
}
n <- 100
sim <- SimData(n)
X_train <- sim$X
y_train <- sim$y
# X_train[1:50,1:10] <- NA
# X_train[1:50, 2] <- NA

sim <- SimData(n)
X_fit <- sim$X
y_true <- sim$y

# In Sample
Trees <- RegTree(y_train,X_train, max_nodes = 31)
tree <- Trees$Tree

y_fit <- FitMat(t(X_train), tree)

ts.plot(cbind(y_train, y_fit), col = c("red", "blue"))

y_true[1]







