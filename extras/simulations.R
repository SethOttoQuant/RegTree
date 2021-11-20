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



run_sim <- function(n){
  sim <- SimData(n)
  X_train <- sim$X
  y_train <- sim$y
  # X_train[1:50,1:10] <- NA
  # X_train[1:50, 2] <- NA
  
  sim <- SimData(n)
  X_fit <- sim$X
  y_true <- sim$y
  
  # --- Regression Nodes -------------------------------------
  # A <- Sys.time()
  Trees <- RegForest(y_train,X_train, max_nodes = 32) # this uses regression at each node
  # Sys.time() - A 

  # Out of Sample
  reg_fit <- FitField(X_fit,Trees) # for RegForest()
  # ts.plot(cbind(fit,y_true), col = c("red", "blue"))

  # --- Standard Nodes  -------------------------------------
  Trees <- RndForest(y_train, X_train, max_nodes = 32)
  std_fit <- StdFitField(X_fit, Trees)
  
  # X_train[1:50,2] <- mean(X_train[51:100,2]) 
  # X_train[1:50, 1:10] <- matrix(1,50,1)%x%t(colMeans(X_train[51:100, 1:10]))
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
  
  res <- data.frame(reg_fit, std_fit,rf_fit,ranger_fit$predictions,y_true)
  names(res) <- c("regression leaves", "standard leaves", "randomForest", "ranger", "true values")
  return(res)
}

A <- Sys.time()
Out <- lapply(rep(500, 10), FUN = run_sim) # 200 observations, 1000 simulations (it's slow) --- 2 hours
B <- Sys.time()
B - A

NewAlgo <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "regression leaves"))
RndFst <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "standard leaves"))
randForest <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "randomForest"))
rang <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "ranger"))
TrueVals <- do.call("cbind", lapply(Out, FUN = get_from_list, what = "true values"))


mean((RndFst - TrueVals)^2)
mean((NewAlgo - TrueVals)^2)
mean((randForest - TrueVals)^2)
mean((rang - TrueVals)^2)



mean((RndFst - randForest)^2)
mean((RndFst - rang)^2)
mean((randForest - rang)^2)


# x = NULL, X, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = ""
pretty_plot(X = res, xlab = "Time", ylab = "Simulated Data", title = "New Algorithm against Existing Libraries")

dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/algo_test.png", width = 1400, height = 900, res = 175)


