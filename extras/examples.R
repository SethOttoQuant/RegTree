library(RegTree)
library(randomForest)
library(ranger)
library(dateutils)

# ------ Simulate Data --------------
SimData <- function(n = 100){
  X <- matrix(rnorm(n*20),n,20)
  b <- c(1,5,3,0,2,-1,2,-3,0,1,1,0,2,0,3,-1,-2,-1,2,4)
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
X_train[1:100,1:10] <- NA
# X_train[1:50, 2] <- NA
sim <- SimData(n)
X_fit <- sim$X
y_true <- sim$y
# fast_cut(X_train[,19], y_train)
# bob <- bestsplit(X_train, y_train, seq(0,NROW(X_train)-1), 10)
# out <- RegTree(y_train, X_train)
# A <- Sys.time()
# out <- RForest(y_train, X_train)
# Sys.time() - A
# fit <- fitfield(X_train, out)
# sum((y_train-fit)^2)

out <- reg_forest(y_train, X_train)

pretty_plot(cbind(out$fit, out$true_vals))

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
  #Trees <- RegForestM(y_train,X_train, depth_range = c(50,100)) # this uses regression at each node
  Trees <- Rnd_Forest(y_train,X_train, depth_range = c(50,100)) # this uses mean(y) at each node
  
  # Out of Sample
  fit <- FitFieldM(X_fit,Trees) # for RegForest()
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
  
  res <- data.frame(fit,rf_fit,ranger_fit$predictions,y_true)
  names(res) <- c("new algo", "randomForest", "ranger", "true values")
  return(res)
}

A <- Sys.time()
Out <- lapply(rep(200, 20), FUN = run_sim) # 200 observations, 100 simulations (it's slow)
B <- Sys.time()
B - A

NewAlgo <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "new algo"))
randForest <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "randomForest"))
rang <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "ranger"))
TrueVals <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "true values"))

mean((NewAlgo - TrueVals)^2, na.rm = TRUE)
mean((randForest - TrueVals)^2)
mean((rang - TrueVals)^2)

mean((NewAlgo - randForest)^2)
mean((NewAlgo - rang)^2)
mean((randForest - rang)^2)






# x = NULL, X, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = ""
# pretty_plot(X = res, xlab = "Time", ylab = "Simulated Data", title = "New Algorithm against Existing Libraries")

# dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/algo_test.png", width = 1400, height = 900, res = 175)


