# ------ Simulate Data --------------
SimData <- function(n = 100){
  X <- matrix(rnorm(n*20),n,20)
  b <- c(1,5,3,0,2,-1,2,-3,0,1,1,0,2,0,3,-1,-2,-1,2,4)
  y <- X%*%b + rnorm(n)
  return(list(X = X,
              b = b,
              y = y))
}

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

sim <- SimData(100)
X_train <- sim$X
y_train <- sim$y

X_train[1:50,2] <- NA

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
rf <- randomForest(x = X_train[,-2], y = c(y_train))
rf_fit <- predict(rf, X_fit[,-2])
ts.plot(cbind(fit,rf_fit,y_true), col = c("red", "blue", "black"))


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