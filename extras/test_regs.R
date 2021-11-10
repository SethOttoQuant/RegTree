x <- 1:40
y <- rep(0,40)
y[1:20] <- x[1:20] + 2*rnorm(20) 
y[21:40] <- -x[21:40] + 4*rnorm(20) 
y[1:10] <- y[1:10] - 10
y[11:20] <- y[11:20] + 10
y[21:30] <- y[21:30] + 10

X <- cbind(rnorm(40), rnorm(40), x, rnorm(40), rnorm(40))
rorder <- select_rnd(40,40) + 1
y <- y[rorder]
X <- X[rorder, ]

oot <- best_split(X,y,1:40,5)
par <- oot[[1]]
idx <- oot[[4]]
er <- oot[[2]]

tst <- y[idx] - par[2] - X[idx,3]*par[3]



cbind(tst, er)


leq <- oot[[4]]
X[leq, 3]

X[oot[[5]],3]



out <- RegTree(y,X,7)
y_fit <- FitMat(t(X), out$Tree)
ts.plot(cbind(y_fit, y), col = c("red", "blue"))


tree <- out$Tree

j <- 11:20

y_fit[j]
y[j]
x[j]

res <- y[1:20] - tree[1,6]
jim <- quickreg(x[1:20], res, 0)
jim[[1]]


FitVec(X[20,], tree)

j <- 20

ft <- rep(tree[1,6],length(j)) # unconditional mean
ft
ft <- ft + tree[2,6] + x[j]*tree[2,7]
ft
ft <- ft + tree[5,6] + x[j]*tree[5,7]
ft

plot(x[j],y[j])
lines(x[j], ft)

ft

res <- y - tree[1,6]
res[]




out <- best_split(X, y, 1:40, 5)
pars <- out[[1]]
pars

ct <- pars[9]
y_fit <- rep(0,40)
y_fit[x<ct] <- pars[2] + pars[3]*x[x<ct]
y_fit[x>ct] <- pars[4] + pars[5]*x[x>ct]
plot(x,y)
lines(x,y_fit,col = "red")

cbind(y[x<ct] - y_fit[x<ct], out[[2]])
cbind(y[x>ct] - y_fit[x>ct], out[[3]])
out[[5]]

tst <- find_cut(x, y, 1:40)
bob <- tst[[1]]
sum(bob[5:7])








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
n <- 200
sim <- SimData(n)
X_train <- sim$X
y_train <- sim$y
X_train[1:100,1:10] <- NA
# X_train[1:50, 2] <- NA

sim <- SimData(n)
X_fit <- sim$X
y_true <- sim$y

# In Sample
# Trees <- RegTree(y_train,X_train, max_nodes = 31)
# tree <- Trees$Tree
# y_fit <- FitMat(t(X_train), tree)
# ts.plot(cbind(y_train, y_fit), col = c("red", "blue"))

Trees <- RegForest(y_train,X_train, max_nodes = 31, draws = 500)
fit <- FitField(X_train,Trees)
ts.plot(cbind(fit,y_train), col = c("red", "blue"))

# out of sample:
fit <- FitField(X_fit,Trees)
ts.plot(cbind(fit,y_true), col = c("red", "blue"))

library(randomForest)
X_train[1:100,1:10] <- 0
rf <- randomForest(x = X_train, y = c(y_train))
rf_fit <- predict(rf, X_fit)
ts.plot(cbind(rf_fit,y_true), col = c("red", "blue"))

mean((fit-y_true)^2)
mean((rf_fit-y_true)^2)













