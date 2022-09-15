library(RegTree)
library(randomForest)
library(ranger)
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

SimDataNonLin <- function(n = 100, has_na=FALSE){
  X <- matrix(rnorm(n*20),n,20)
  b1 <- c( 5, 1,-5,-1, 2,-1, 0,-3, 2, 1, 1, 0, 0, 0, 2,-1, 2,-1, 0, 1)
  b2 <- c(-5,-2, 0, 2,-1,-2, 4, 3, 1,-1, 1, 0,-2, 0, 1,-1,-2, 1, 2, 0)
  b3 <- c( 1, 0, 5, 0,-2, 3,-1, 0,-1, 1, 0,-1, 0, 1,-2, 0, 0, 0, 1,-1)
  z <- rnorm(n)
  # zx <- runif(n)
  # X[zx<1/2,11:20] <- X[zx<1/2,11:20] + 2
  # X[zx>1/2,1:10] <- X[zx>1/2,1:10] - 5
  y <- X%*%b1
  y[z < -1] <- X[z < -1, ]%*%b2 + 1
  if(has_na){
    X[z < -2, 1:10] <- NA
  }
  y[z > 1] <- X[z > 1, ]%*%b3 - 1 
  if(has_na){
    X[z>3,1:10] <- NA
  }
  # X <- exp(X)
  # X[z<1/3,11:20] <- X[z<1/3,11:20] - 1
  # y <- exp(y/10) + rnorm(n)/4
  y <- y + rnorm(n)
  X <- cbind(c(z) + rnorm(n)/4, X)
  return(list(X = X,
              y = y))
}

# SimData3 <- function(n = 100){
#   X <- matrix(rnorm(n*10),n,10)
#   b3 <- c( 1, 0, 5, 0,-2, 3,-1, 0,-1, 1)
#   y <- X%*%b3
#   y <- y + rnorm(n)
#   # X <- cbind(X,c(z))
#   return(list(X = X,
#               y = y))
# }

n <- 600
# Run a single simulation
sim <- SimDataNonLin(n)
X <- sim$X
y <- sim$y
X[1:100,1:3] <- NA

A <- Sys.time()
out_std <- RndForest(y,X, max_obs = 10, min_obs = 3) # this uses mean(y) at each node
# tst <- reg_forest(y, X, max_obs = 15, min_obs = 3, weight_by_mse = FALSE, type = "standard")
Sys.time() - A 

sim <- SimDataNonLin(200)
X_fit <- sim$X
y_true <- sim$y

fit_std <- StdFitField(X_fit,out_std$Trees, weight_nodes = TRUE)

pretty_plot(cbind(y_true, fit_std[[1]]))
mean((y_true - fit_std[[1]])^2)

# Tree <- out_std$Trees[[9]]
# Tree <- cbind(Tree, Tree[,7]/(Tree[,8]-1))
# sue <- fitmat(t(X_fit), out_std$Trees[[9]], weight = TRUE)

# x <- c(3,2,5,4,6,3,2,3,4,1)
# y <- c(1,3,2,4,5,6,2,7, 5)
# 
# m_1 <- mean(x)
# m_2 <- mean(y)
# s_1 <- sum((x - m_1)^2)
# s_2 <- sum((y - m_2)^2)
# n_1 <- length(x)
# n_2 <- length(y)
# 
# welch_t(m_1, m_2, s_1, s_2, n_1, n_2)
# out <- t.test(y, x, "less")
# pt(out$statistic, out$parameter, lower.tail = FALSE)
# pt(1, 3)
# pnorm(0)
# t_tst(1,3)


# out <- StdFitFieldWeight(X_fit, tst$Trees)
n <- 200

run_sim <- function(n){
  sim <- SimDataNonLin(n)
  # sim <- SimData(n)
  X_train <- sim$X
  y_train <- sim$y
  # X_train[1:100,1:10] <- NA
  # X_train[1:50, 2] <- NA
  
  sim <- SimDataNonLin(n, has_na = FALSE)
  # sim <- SimData(n)
  X_fit <- sim$X
  y_true <- sim$y
  
  # In Sample
  out_std <- RndForest(y_train,X_train, max_obs = 15, min_obs = 5) # this uses mean(y) at each node
  out_reg <- RegForest(y_train,X_train, min_obs = 15, max_nodes = 100) 
  
  # Out of Sample
  fit_reg <- FitField(X_fit,out_reg$Trees)[[1]] # for RegForest()
  fit_std <- StdFitField(X_fit,out_std$Trees, weight_nodes = FALSE)[[1]]
  # fit2 <- StdFitField(X_fit,Trees2) # for RndForest()
  # ts.plot(cbind(fit_std,y_true), col = c("red", "blue"))
  # isna <- is.na(X_train[,1])
  # X_train[1:100, 1:10] <- matrix(1,sum(isna),1)%x%t(colMeans(X_train[, 1:10, drop=FALSE], na.rm=TRUE))
  # X_train <- X_train[51:100,]
  # y_train <- y_train[51:100]
  # X_train <- X_train[ , -seq(10)]
  # X_fit <- X_fit[,-seq(10)]
  
  # vs randomForest package
  rf <- randomForest(x = X_train, y = c(y_train))
  rf_fit <- predict(rf, X_fit)
  
  df <- data.frame(y_train, X_train)
  rf2 <- ranger(y_train ~ ., data = df)
  ranger_fit <- predict(rf2, data = data.frame(X_fit))
  
  res <- data.frame(fit_reg, fit_std, rf_fit, ranger_fit$predictions, y_true)
  names(res) <- c("regression", "standard", "randomForest", "ranger", "true values")
  return(res)
}

A <- Sys.time()
Out <- lapply(rep(200, 1000), FUN = run_sim) # 200 observations, xxx simulations (it's slow)
B <- Sys.time()
B - A

Reg <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "regression"))
Std <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "standard"))
randForest <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "randomForest"))
rang <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "ranger"))
TrueVals <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "true values"))

mean((Reg - TrueVals)^2, na.rm = TRUE)
mean((Std - TrueVals)^2, na.rm = TRUE)
mean((randForest - TrueVals)^2)
mean((rang - TrueVals)^2)

j <- 2
pretty_plot(cbind(TrueVals[,j], Std[,j], rang[,j]))


