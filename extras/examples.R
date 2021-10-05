library(RegTree)
library(randomForest)
library(ranger)

# ------ Simulate Data --------------
SimData <- function(n = 100){
  X <- matrix(rnorm(n*20),n,20)
  b <- c(1,10,3,0,2,-1,2,-3,0,1,1,0,2,0,3,-1,-2,-1,2,4)
  y <- X%*%b + rnorm(n)
  return(list(X = X,
              b = b,
              y = y))
}

sim <- SimData(100)
X_train <- sim$X
y_train <- sim$y

sim <- SimData(100)
X_fit <- sim$X
y_true <- sim$y

# In Sample
Trees <- RegForest(y_train,X_train, max_nodes = 32)
# fit <- FitField(X_train,Trees)
# ts.plot(cbind(fit,y_train), col = c("red", "blue"))

# Out of Sample
fit <- FitField(X_fit,Trees)
# ts.plot(cbind(fit,y_true), col = c("red", "blue"))

# vs randomForest package
rf <- randomForest(x = X_train, y = c(y_train))
rf_fit <- predict(rf, X_fit)

df <- data.frame(y_train, X_train)
rf2 <- ranger(y_train ~ ., data = df)
ranger_fit <- predict(rf2, data = data.frame(X_fit))

res <- data.frame(fit,rf_fit,ranger_fit$predictions,y_true)
names(res) <- c("new algo", "randomForest", "ranger", "true values")

# x = NULL, X, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = ""
pretty_plot(X = res, xlab = "Time", ylab = "Simulated Data", title = "New Algorithm against Existing Libraries")

dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/algo_test.png", width = 1400, height = 900, res = 175)

X_train[1:50,2] <- NA

# In Sample
Trees <- RegForest(y_train,X_train, max_nodes = 32)
# fit <- FitField(X_train,Trees)
# ts.plot(cbind(fit,y_train), col = c("red", "blue"))

# Out of Sample
fit <- FitField(X_fit,Trees)
# ts.plot(cbind(fit,y_true), col = c("red", "blue"))

# vs randomForest package
rf <- randomForest(x = X_train[,-2], y = c(y_train))
rf_fit <- predict(rf, X_fit[,-2])

df <- data.frame(y_train, X_train[,-2])
rf2 <- ranger(y_train ~ ., data = df)
ranger_fit <- predict(rf2, data = data.frame(X_fit[,-2]))

res <- data.frame(fit,rf_fit,ranger_fit$predictions,y_true)
names(res) <- c("new algo", "randomForest", "ranger", "true values")

# x = NULL, X, lwd = 2, xlab = "Date", ylab = "", legend_pos = "bottomleft", title = ""
pretty_plot(X = res, xlab = "Time", ylab = "Simulated Data", title = "New Algorithm against Existing Libraries")

dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/algo_test_na.png", width = 1400, height = 900, res = 175)










# --- Old Stuff -------------------------

library(devtools)
load_all("~/RegTree")

sim <- SimData()
X <- sim$X
y <- sim$y
select_rnd(10,5)



tmp <- best_split(X,y,10)

out <- RegTree(y,X, 64, .02)

X[1:50,1] <- NA
X[51:100,2] <- NA

out <- RegTree(y,X, 64, .02)
Fit <- FitMat(t(X),Tree = out)
ts.plot(cbind(Fit,y), col = c("red", "blue"))




getTree(rf, 40)




FitMat(t(X),Trees[[2]])


Trees[[2]]



500*(B-A)


fit <- FitVec(X[100, ], out)

A <- Sys.time()
for(j in 1:1000){
  Fit <- apply(X, 1, FitVec, Tree = out)
}
B <- Sys.time()
B-A

A <- Sys.time()
for(j in 1:1000){
  Fit2 <- FitMat(t(X),Tree = out)
}
B <- Sys.time()
B-A

FitVec(X[3,], out)

bob <- cbind(Fit,Fit2)

ts.plot(cbind(Fit,y), col = c("red", "blue"))


cdn <- node_conditions(out, 6)

idx <- X[, 18] > -0.2964146 & X[ ,1] < 0.2870587 & X[ ,20] > -0.6341549

mean(y[idx])

idx <- X[,7] < 1.3185637
a <- sum((y[idx] - mean(y[idx]))^2)
b <- sum((y[!idx] - mean(y[!idx]))^2)


mean(y)

out