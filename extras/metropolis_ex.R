library(RegTree)
library(dateutils)
library(devtools)
# load_all()

# ------ Simulate Data --------------
SimData <- function(n = 100){
  X <- matrix(rnorm(n*20),n,20)
  b <- c(1,5,3,0,2,-1,2,-3,0,1,1,0,2,0,3,-1,-2,-1,2,4)
  y <- X%*%b + rnorm(n)
  return(list(X = X,
              b = b,
              y = y))
}

x <- c(0,1,2,3,2,1,4,6,7,8,9,10)
y <- c(1,2,1,2,3,4,2,6,5,6,5,5)
select_cut(x,y,0,4,1)


n <- 200
# Run a single simulation
train <- SimData(n)
X <- train$X
y <- train$y
# tst <- Rsimtree(y, X)
# tst[[1]]
# tst[[2]][[20]]
test <- SimData(n)
# X[1:100,1:3] <- NA

out <- reg_forest(y, X, max_nodes = 41, draws = 5000, type = "standard")

tree <- out$Trees[[1]]

fit <- StdFitField(test$X, out$Trees)
fit <- fit[[1]]
pretty_plot(cbind(test$y, fit))



 # ------------------------------



out <- simforest(y, X, max_nodes = 21, burn = 200, reps=5000)

tree <- out$Trees[[2]]






tree <- Trees[[1]]

v <- c(1,2,3,5,6,8,9,12)
find_unique(v,12)

