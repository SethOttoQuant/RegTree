rm(list=ls()) # clear workspace

# Packages -----------------------------

library(randomForest)
library(ranger)

library(RegTree)
library(dateutils)

library(tidyverse)
library(tsbox)
library(forecast)
library(tsDyn)

library(future.apply)
plan(multisession)
library(parallel)
library(foreach)
library(doParallel)
n_cores <- detectCores()
cl <- makeCluster(n_cores)
registerDoParallel(cores = detectCores())

library(MonteCarlo)
library(Metrics)

# Medeiros Replication
library(HDeconometrics)
library(glmnet)
library(fbi)
library(TTR)


# 1) Simple linear Model ----
Model.linear <- function(n_t = 100, n_x = 10, err_sd = 1){
  set.seed(1234)
  X <- matrix(rnorm(n_t*n_x, mean = 0, sd = 1), n_t, n_x)
  # b <- sample(-1:1, n_x, replace=T)
  b <- runif(n_x, min = -1, max = 1)
  y <- X%*%b + rnorm(n_t, mean = 0, sd = err_sd)
  return(list(X = ts(X),
              b = b,
              y = ts(y)))
}

# parameter grids
n_t <- c(100, 200, 400, 800)
n_x <- c(10, 25, 100, 200)

# Benchmark is n_t[3]; n_x[2]

err_sd <- 0.95

reps <- 50 # monte carlo replications
train_frac <- 0.75 # split training sample and test sample
contains_na <- "FALSE"
num_na <- "one" # options: "one", "frac"
na_xfrac <- 0.2
na_tfrac <- 0.5
replace_by <- "drop" # options: "mean", "drop"

# # Grid search for parameter detection
# 
# sd_err <- seq(from = 0.01, to = 50, by = 0.01)
# 
# for(i in 1:length(sd_err)){
#   out <- Model.linear(n_t[2], n_x[2], err_sd = sd_err[i])
#   fit1 <- lm(out$y ~ out$X)
#   corr <- cor(fitted(fit1), out$y )^2
#   if (round(corr, digits = 2) == 0.80){
#     print(sd_err[i])
#     break
#   }
# }
# 
# # Standard Errors
# # T     |   X     |  STD 
# # -----------------------
# # 400   |   25    |  1.51 
# # -----------------------
# # 100   |   25    |  1.97 
# # 200   |   25    |  1.65
# # 800   |   25    |  1.58 
# # -----------------------
# # 400   |   10    |  0.91 
# # 400   |   100   |  3.48 
# # 400   |   200   |  8.35

# simulation function  
lin_mod_boot <- function(n_t, 
                         n_x, 
                         err_sd, 
                         train_frac, 
                         contains_na,
                         num_na,
                         na_xfrac,
                         na_tfrac,
                         replace_by
                         ) {
  
  if(n_t == 400 & n_x == 25){
    err_sd = 1.51} else if (n_t == 100 & n_x == 25){
      err_sd = 1.97} else if (n_t == 200 & n_x == 25){
        err_sd = 1.65} else if (n_t == 800 & n_x == 25){ 
          err_sd = 1.58} else if (n_t == 400 & n_x == 10){
            err_sd = 0.91} else if (n_t == 400 & n_x == 100){
              err_sd = 3.48} else if (n_t == 400 & n_x == 200){
                err_sd = 8.35} 
  
  out.sim <- Model.linear(n_t, n_x, err_sd = err_sd)
  
  train.idx <- sample(nrow(out.sim$y), train_frac * nrow(out.sim$y))
  x.train <- out.sim$X[train.idx, ]
  y.train <- out.sim$y[train.idx, ]
  
  if(isTRUE(contains_na)){
    if (num_na == "one"){
      x.train[1:(na_tfrac*n_t),1] <- NA
    } else {
      x.train[1:(na_tfrac*n_t), 1:(na_xfrac*n_x)] <- NA
    }
  }
  
  x.test <- out.sim$X[-train.idx, ]
  y.test <- out.sim$y[-train.idx, ]
  
  RMSE <- list()
  # # Linear Prediction ARIMA(0,0,0) -----
  # fit <- lm(y.train ~ x.train)
  # 
  # fcasts <- matrix(NA, length(y.test),1)
  # for (i in 1:length(y.test)) { # start rolling forecast
  #   # start from 1997, every time one more year included
  #   if(i == 1){
  #     fcasts[[i]] <- forecast(fit, xreg = x.train, h = 1)$mean[1]
  #   } else {
  #     fcasts[[i]] <- forecast(fit, xreg = rbind(x.test, x.test[1:i,]))$mean[i]
  #   }
  # }

  ### metric ####
  # RMSE$linear <- rmse(y.test, fcasts) # return RMSE

  # RF own ------------------------------
  rf_own <- RndForest(y.train, 
                      x.train, 
                      max_obs = 3, # max size of terminal nodes
                      min_obs = 1, # min size of terminal nodes 
                      max_nodes = 2000 # max number of terminal nodes trees in the forest can have
  )
  # Make predictions based on trained forest
  fcasts <- matrix(NA, length(y.test),1)
  for (i in 1:length(y.test)) { # start rolling forecast
    if(i == 1){
      pred.rf <- StdFitField(as.matrix(rbind(x.train, x.test[1,])), rf_own$Trees)[[1]]
      fcasts[[i]] <- tail(pred.rf,1)
    } else {
      # pred.rf <- StdFitField(as.matrix(rbind(x.train, x.test[1:i,])), rf_own$Trees)[[1]]
      pred.rf <- StdFitField(as.matrix(x.test[1:i,]), rf_own$Trees)[[1]]
      fcasts[[i]] <- tail(pred.rf,1)
    }
  }
  RMSE$rf <- rmse(y.test, fcasts) # return RMSE  
   
  # RFRN own ----------------------------
  rfrn_own <- RegForest(y.train, 
                        x.train, 
                        min_obs = 10, # min size of terminal nodes
                        max_nodes = 2000 # max number of terminal nodes trees in the forest can have.
  )
  
  # Make predictions based on trained forest
  fcasts <- matrix(NA, length(y.test),1)
  for (i in 1:length(y.test)) { # start rolling forecast
    if(i == 1){
      pred.rfrn_own <- FitField(as.matrix(rbind(x.train, x.test[1,])), rfrn_own$Trees)[[1]]
      fcasts[[i]] <- tail(pred.rf,1)
    } else {
      # pred.rfrn_own <- FitField(as.matrix(rbind(x.train, x.test[1:i,])), rfrn_own$Trees)[[1]]
      pred.rfrn_own <- FitField(as.matrix(x.test[1:i,]), rfrn_own$Trees)[[1]]
      fcasts[[i]] <- tail(pred.rf,1)
    }
  }
  RMSE$rfrn <- rmse(y.test, fcasts) # return RMSE 
  
  # Handling Missings ----------------
  
  if(isTRUE(contains_na)){
    if (num_na == "one"){
      if (replace_by == "mean"){
        x.train[1:(na_tfrac*n_t),1] <- matrix(1,(na_tfrac*n_t),1)%x%t(colMeans(x.train[, 1, drop=FALSE], na.rm=TRUE))
      } else {
        y.train <- y.train[!is.na(x.train[,1])]
        x.train <- na.omit(x.train)
      }
    } else {
      if (replace_by == "mean"){
        x.train[1:(na_tfrac*n_t),1:(na_xfrac*n_x)] <- matrix(1,(na_tfrac*n_t),(na_xfrac*n_x))%x%t(colMeans(x.train[, (na_xfrac*n_x), drop=FALSE], na.rm=TRUE))
      } else {
        y.train <- y.train[!is.na(x.train[,1])]
        x.train <- na.omit(x.train)
      }
    }
  }
    
  # Ranger ------------------------------
  df <- data.frame(y.train, x.train)
  # Train the forest
  rg <- ranger(y.train ~ ., 
               data = df,
               min.node.size = 10)
  
  # Make predictions based on trained forest
  fcasts <- matrix(NA, length(y.test),1)
  for (i in 1:length(y.test)) { # start rolling forecast
    if(i == 1){
      fcasts[[i]] <- predict(rg, data = data.frame(rbind(x.train, x.test[1,])))$predictions[2]
    } else {
      # fcasts[[i]] <- predict(rg, data = data.frame(rbind(x.train, x.test[1:i,])))$predictions[i]
      fcasts[[i]] <- predict(rg, data = data.frame(x.test[1:i,]))$predictions[i]
    }
  }
  RMSE$ranger <- rmse(y.test, fcasts) # return RMSE
  
  # RANDOM FOREST -------------------------
  rf <- randomForest(x = x.train, 
                     y = c(y.train),
                     nodesize = 10
  )
  # Make predictions based on trained forest
  fcasts <- matrix(NA, length(y.test),1)
  for (i in 1:length(y.test)) { # start rolling forecast
    if(i == 1){
      pred.rf <- predict(rf, as.matrix(rbind(x.train, x.test[1,])))
      fcasts[[i]] <- tail(pred.rf,1)
    } else {
      # pred.rf <- predict(rf, as.matrix(rbind(x.train, x.test[1:i,])))
      pred.rf <- predict(rf, as.matrix(x.test[1:i,]))
      fcasts[[i]] <- tail(pred.rf,1)
    }
  }
  RMSE$randfor <- rmse(y.test, fcasts) # return RMSE

  # MEDEIROS ET AL -------------------------
  
  x0 = as.matrix(rbind(x.train))
  
  factors = princomp(scale(x0))$scores[,1:4]
  x=cbind(x0,factors)
  X=embed(as.matrix(x),4)
  
  Xin=X[-c((nrow(X)):nrow(X)),]
  yin=tail(y.train,nrow(Xin))
  
  rf_medeiros = randomForest::randomForest(Xin, yin, importance = TRUE)
  
  fcasts <- matrix(NA, length(y.test),1)
  for (i in 1:length(y.test)) { # start rolling forecast
    if(i == 1){
      
      x0 = as.matrix(rbind(x.train, x.test[1,]))
      
      factors = princomp(scale(x0))$scores[,1:4]
      x=cbind(x0,factors)
      X=embed(as.matrix(x),4)
      
      Xout=X[nrow(X),]
      Xout=t(as.vector(Xout))
      
      fcasts[[i]] <- predict(rf_medeiros,Xout)
    } else {
      
      x0 = as.matrix(rbind(x.train, x.test[1:i,]))
      
      factors = princomp(scale(x0))$scores[,1:4]
      x=cbind(x0,factors)

      X=embed(as.matrix(x),4)
      
      Xout=X[nrow(X),]
      Xout=t(as.vector(Xout))
      
      fcasts[[i]] <- predict(rf_medeiros,Xout)
    }
  }
  
  RMSE$medeiros <- rmse(y.test, fcasts) # return RMSE
  
  return(
    list("Ranger" = RMSE$ranger,
         "Randfor" = RMSE$randfor,
         "RF" = RMSE$rf, 
         "RFRN" = RMSE$rfrn,
         "Medeiros" = RMSE$medeiros)
  )
}

param_list = list("n_t" = n_t[3],
                  "n_x" = n_x,
                  "err_sd" = err_sd, 
                  "train_frac" = train_frac,
                  "contains_na" = contains_na,
                  "num_na" = num_na, 
                  "na_xfrac" = na_xfrac,
                  "na_tfrac" = na_tfrac,
                  "replace_by" = replace_by)

set.seed(1234)

MC_result.linear <- MonteCarlo(func = lin_mod_boot, 
                        nrep = reps,
                        ncpus = parallel::detectCores() - 1,
                        param_list = param_list,
                        export_also = list(
                          "packages" = c("forecast", "Metrics")
                        ),
                        time_n_test = TRUE,
                        raw = T)

Frame <- MakeFrame(MC_result.linear)
Frame 

colMeans(Frame[c("Linear", "Ranger")])
