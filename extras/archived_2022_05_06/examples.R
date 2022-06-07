# library(RegTree)
# library(randomForest)
# library(ranger)
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

SimDataNonLin <- function(n = 100){
  X <- matrix(rnorm(n*10),n,10)
  b1 <- c(1,5,3,0,2,-1,0,-3,0,0)
  b2 <- c(-2,0,0,2,-1,0,4,-3,1,-1)
  b3 <- c(0,0,1,0,0,-3,3,0,0,1)
  z <- rnorm(n)
  y <- X%*%b1
  y[z<(-.5)] <- X[z<(-.5), ]%*%b2
  y[z>.5] <- X[z>.5, ]%*%b3
  y <- y + rnorm(n)
  X <- cbind(X,z)
  return(list(X = X,
              y = y))
}

n <- 200
# Run a single simulation
sim <- SimDataNonLin(n)
X <- sim$X
y <- sim$y
X[1:100,1:3] <- NA
sim <- SimDataNonLin(n)
X_fit <- sim$X
y_true <- sim$y

A <- Sys.time()
tst <- reg_forest(y, X, weight_by_mse = TRUE, type = "alt", geom_par = .5, weight_pow = 8)

out <- StdFitFieldWeight(X_fit, tst$Trees, tst$)

ts.plot(cbind(tst$true_vals, tst$out_of_sample), col = c("blue", "red"))
tst$mse

tst$MSE

rand_geom(.5)


bob <- Reg_Forest(y, X)

oob <- rowMeans(bob$OOB, na.rm = TRUE)

tmp <- FitMat(t(X), bob$Trees[[2]])
tmp[[1]]



jim <- FitField(X, bob$Trees)
jim[[1]]

sue <- FitVec(X[2,], bob$Trees[[2]])
sue[[1]]

# fast_cut(X_train[,19], y_train)
# bob <- bestsplit(X_train, y_train, seq(0,NROW(X_train)-1), 10)
# out <- RegTree(y_train, X_train)
# A <- Sys.time()
# out <- RForest(y_train, X_train)
# Sys.time() - A
# fit <- fitfield(X_train, out)
# sum((y_train-fit)^2)

to_keep <- select_rnd(n, 120);
oout <- regtree(y_train, X_train, to_keep[[1]])

jim <- fitvec(X_fit[1,], oout)
carl <- jim[[2]]
sum(carl) + oout[1,6]
jim[[1]]



jane <- fitmat(t(X_fit), oout)
jen <- jane[[2]]

colSums(jen) + oout[1,6]
c(jane[[1]] )

tst <- rforest(y_train, X_train)
bob <- tst$FC[,,1]


out <- Reg_Forest(y_train, X_train)
FC <- out$FC

FCm <- apply(FC, c(1,2), FUN = mean, na.rm=TRUE)
out = FitMat(t(X_train), Tree)
ft = out[[1]]
fc = t(out[[2]])

tmp <- FitField(X_fit, out$Trees)




Tree = Reg_Tree(y_train, X_train, seq(0,199), 31)



mean(y_train[to_keep[[1]]+1])



rs = rowSums(fc)

j = 1
rs[j]
ft[j]
# tst = FitMat(X_fit, Tree)











FitVec(X_train[20,],Tree)

out <- reg_forest(y_train, X_train, regression = TRUE)

pretty_plot(cbind(y_train, out$out_of_sample))

out <- reg_forest(y_train, X_train)
pretty_plot(cbind(out$true_vals, out$out_of_sample))

# A <- Sys.time()
# Trees <- rand_rand(X_train, c(50,100), 10000)
# Sys.time() - A
# 
# out <- rnd_fit(X_fit, X_train, y_train, Trees, 1)
# 
# ts.plot(cbind(out[[1]], y_true), col = c("red", "blue"))
# 
# bob <- fitmat(t(X_fit), X_train, y_train, Trees[[2]])
# 
# jim <- fitvec(X_fit[33,], X_train, y_train, Trees[[2]])
# 
# 
# bob <- do.call("c", sapply(Trees, function(j) j[,6]))
# any(is.na(bob))
# 
# B <- Sys.time()
# Trees <- Rnd_Forest(y_train,X_train, 30)
# Sys.time()-B

# B <- Sys.time()
# Trees <- Reg_Forest(y_train,X_train, 30)
# Sys.time()-B

run_sim <- function(n){
  sim <- SimData(n)
  X_train <- sim$X
  y_train <- sim$y
  X_train[1:100,1:10] <- NA
  # X_train[1:50, 2] <- NA
  
  
  sim <- SimData(n)
  X_fit <- sim$X
  y_true <- sim$y
  
  # In Sample
  out_std <- RndForest(y_train,X_train)
  out_reg <- RegForest(y_train,X_train) # this uses mean(y) at each node
  
  # Out of Sample
  fit_reg <- FitField(X_fit,out_reg$Trees) # for RegForest()
  fit_std <- StdFitField(X_fit,out_std$Trees)
  # fit2 <- StdFitField(X_fit,Trees2) # for RndForest()
  # ts.plot(cbind(fit,y_true), col = c("red", "blue"))
  
  # X_train[1:50,2] <- mean(X_train[51:100,2]) 
  X_train[1:100, 1:10] <- matrix(1,100,1)%x%t(colMeans(X_train[101:200, 1:10]))
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
Out <- lapply(rep(300, 20), FUN = run_sim) # 200 observations, 100 simulations (it's slow)
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


