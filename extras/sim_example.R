# example using predict for ranger:

library(ranger)
library(randomForest)
library(RegTree)
set.seed(1234)

Model.linear <- function(n_t = 100, n_x = 10, err_sd = 1, x_sd = 1, b = NULL){
  if(!is.null(b)){
    if(length(b)!=n_x){
      stop("Dimensions of 'b_in' do not agree with 'n_x'")
    }
  }else{
    b <- runif(n_x, min = -1, max = 1)
  }
  X <- matrix(rnorm(n_t*n_x, mean = 0, sd = x_sd), n_t, n_x)
  y <-ts(X%*%b + rnorm(n_t, mean = 0, sd = err_sd))
  colnames(y) <- "y"
  return(list(X = ts(X),
              b = b,
              y = y))
}

# train on 400 observations and test on 100 observations
# testing the predict() function for ranger
n_x <- 10
n_t <- 400
b <- runif(n_x, min = -1, max = 1) # fix b so it doesn't change between draws
training <- Model.linear(n_t = n_t, n_x = n_x, err_sd = 1, x_sd = 1, b = b)
testing <- Model.linear(n_t = 100, n_x = n_x, err_sd = 1, x_sd = 1, b = b)
df <- data.frame(training$y, training$X)
rg <- ranger(y ~ ., data = df, min.node.size = 10)
outA <- predict(rg, data = data.frame(testing$X))
pred_A <- outA$predictions
outB <- predict(rg, data.frame(rbind(training$X, testing$X)))
pred_B <- tail(outB$predictions, 100)
max(pred_A - pred_B)

# Comparing forecasts when the variance of X increases in the testing data

# train on 400 observations and test on 100 observations
# testing the predict() function for ranger
n_x <- 5
n_t <- 400
b <- runif(n_x, min = -1, max = 1) # fix b so it doesn't change between draws
training <- Model.linear(n_t = n_t, n_x = n_x, err_sd = 1, x_sd = 1, b = b)
testing <- Model.linear(n_t = 100, n_x = n_x, err_sd = 1, x_sd = 2, b = b)
df <- data.frame(training$y, training$X)
rg <- ranger(y ~ ., data = df, min.node.size = 10)
pred_rg <- predict(rg, data = data.frame(testing$X))

rf <- reg_forest(training$y, training$X, type = "regression")
pred_rf <- FitField(testing$X, rf$Trees)

plot_data <- data.frame("True" = testing$y, "ranger" = pred_rg$predictions, "RF_RN" = pred_rf[[1]])
pretty_plot(plot_data)

