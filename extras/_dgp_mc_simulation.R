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
library(MSwM)

# DGPs ---------------------

# 1) Simple linear Model ----
Model.linear <- function(n_t = 100, n_x = 10, err_sd = 1){
  set.seed(1234)
  X <- matrix(rnorm(n_t*n_x, mean = 0, sd = 1), n_t, n_x)
  # b <- sample(-1:1, n_x, replace=T)
  b <- runif(n_x, min = -1, max = 1) # Did we want to keep these numbers the same between runs? 
  y <- X%*%b + rnorm(n_t, mean = 0, sd = err_sd)
  return(list(X = ts(X),
              b = b,
              y = ts(y)))
}

# 2) ARIMAX-Model ---- 
Model.arimax <- function(n_t = 100, n_x = 10, rho = 0.9, err_sd_x = 1, err_sd_y = 1){
  set.seed(1234)
  X <- replicate(n_x, arima.sim(n = n_t, model = list(ar = rho), sd = err_sd_x))
  # X <- matrix(rnorm(n_t*n_x, mean = 0, sd = 1), n_t, n_x)
  # b <- sample(-5:5, n_x, replace=T)
  b <- runif(n_x, min = -1, max = 1)
  y_ar <- arima.sim(list(ar = rho), n = n_t, sd = err_sd_y)
  y <- y_ar + X%*%b 
  return(list(X = ts(X),
              b = b,
              y = ts(y)))
}

# 3) SETAR-Model (self-exciting threshold AR) ----
Model.setar <- function(n_t = 100,        # T
                        n_x = 10,
                        # err_sd = 0.6,
                        alpha = 0.4,  # AR(1) of regime 1
                        gamma = 0.8,    # AR(1) of regime 2
                        th_r = -1.0,          # threshold to be in regime 1
                        sigma = 1,  # Variance of error 
                        delta = 3,  # Size of error
                        burnin = 100){
  
  set.seed(1234)
  X <- matrix(rnorm((n_t+burnin)*n_x, mean = 0, sd = 1), (n_t+burnin), n_x)
  # b <- sample(-1:1, n_x, replace=T)
  b1 <- runif(n_x, min = -1, max = 1)
  b2 <- runif(n_x, min = -3, max = 3)
  # y <- X%*%b + rnorm((n_t+burnin), mean = 0, sd = err_sd)
  
  
  # Generate noise
  e <- rnorm(n_t+burnin, 0, sigma)
  # Create space for y
  y <- numeric(n_t+burnin)
  xb1 <- X%*%b1
  xb2 <- X%*%b2
  # Generate time series
  for(i in 2:(n_t+burnin))
  {
    if(y[i-1] <= th_r)
      y[i] <- alpha*y[i-1] + xb1[i] + e[i]
    else
      y[i] <- gamma*y[i-1] + xb2[i] + delta*e[i]
  }
  # Throw away first burnin values
  y <- ts(y[-(1:burnin)])
  # Return result
  return(list(X = ts(X),
              # b = b,
              y = ts(y)))
}

# 4) Markov Switching Process ----
Model.MS <- function(n_t = 100, 
                     n_x = 10,
                     c0_1 = 0,
                     c0_2 = 0,
                     rho = 0.9,
                     sigma_1 = 1,
                     sigma_2 = 0.5,
                     p11 = 0.95,
                     p22 = 0.85){
  
  P <- matrix(c(p11,1-p22,1-p11,p22),2,2)
  s0 <- 1 
  st <- function(i) sample(1:2,1,prob = P[i,])
  
  s <- st(s0)
  for(t in 2:n_t) {
    s <- c(s,st(s[t-1]))
  }
  
  set.seed(1234)
  
  y <- rep(NA, n_t)
  X <- matrix(rnorm((n_t * n_x), mean = 0, sd = 1), n_t, n_x)
  b <- runif(n_x, min = -1, max = 1)
  
  err_1 <- rnorm(n_t, 0, sigma_1)
  err_2 <- rnorm(n_t, 0, sigma_2)
  
  xb <- X%*%b
  y[1] <- c0_1 + xb[1] + err_1[1]
  
  for (t in 2:n_t){
    if(s[t]==2){
      y[t] <- c0_1 + xb[t] + err_1[t]
    } else if (s[t]==1){
      y[t] <- c0_2 + rho*y[t-1] + err_2[t]  
    }
  }
  return(list(X = ts(X),
              b = b,
              y = ts(y)))
}

# 5 Nonlinear Model ----

Model.nonlinear <- function(n_t = 100,        # T
                            n_x = 10,
                            sigma = 0.6,
                            rho = 0.8, 
                            burnin = 100){
  set.seed(1234)
  X <- matrix(rnorm((n_t+burnin)*n_x, mean = 0, sd = 1), (n_t+burnin), n_x)
  b <- runif(n_x, min = -1, max = 1)
  
  # Generate noise
  e <- rnorm((n_t+burnin), 0, sigma)
  # Create space for y
  y <- numeric((n_t+burnin))
  y[1] <- mean(X%*%b)
  fx <- numeric((n_t+burnin))
  # Generate time series
  for(i in 2:(n_t+burnin)){
    fx[i] <- exp(X[i,])%*%sin(X[i,])
    y[i] <- rho*y[i-1] + fx[i] + e[i]
  }
  # Throw away first burnin values
  
  return(list(X = ts(X[-(1:burnin)]),
              b = b,
              y = ts(y[-(1:burnin)])))
}


# TEST SIMULATIONS -------------------------------------------------

n_t <- c(100, 200, 400, 800)
n_x <- c(10, 25, 100, 200)

# Linear Model
test.lin_5 <- Model.linear(n_t[2], n_x[1], err_sd = 0.95)
test.lin_20 <- Model.linear(n_t[2], n_x[2], err_sd = 1.65)
test.lin_100 <- Model.linear(n_t[2], n_x[3], err_sd = 5.2)

ts_c(`5 Xs` = (test.lin_5$y), 
     `20 Xs` = (test.lin_20$y), 
     `100 Xs` = (test.lin_100$y)) %>% 
  ts_plot()


fit1 <- lm(test.lin_5$y ~ test.lin_5$X)
fit2 <- lm(test.lin_20$y ~ test.lin_20$X)
fit3 <- lm(test.lin_100$y ~ test.lin_100$X)
cor(fitted(fit1),test.lin_5$y )^2
cor(fitted(fit2),test.lin_20$y )^2
cor(fitted(fit3),test.lin_100$y )^2

# ARIMAX Model
test.ar_5 <- Model.arimax(n_t[2], n_x[1], rho = 0.8, err_sd_y = 3, err_sd_x = 1.6)
test.ar_20 <- Model.arimax(n_t[2], n_x[2], rho = 0.8, err_sd_y = 3.8, err_sd_x = 1.1)
test.ar_100 <- Model.arimax(n_t[2], n_x[3], rho = 0.8, err_sd_y = 9, err_sd_x = 0.1)

ts_c(`5 Xs` = (test.ar_5$y), 
     `20 Xs` = (test.ar_20$y), 
     `100 Xs` = (test.ar_100$y)) %>% 
  ts_plot()

ts_c(`5 Xs` = (test.lin_5$y), 
     `5 Xs AR` = (test.ar_5$y)) %>% 
  ts_plot()

ts_c(`20 Xs` = (test.lin_20$y), 
     `20 Xs AR` = (test.ar_20$y)) %>% 
  ts_plot()

ts_c(`100 Xs` = (test.lin_100$y), 
     `100 Xs AR` = (test.ar_100$y)) %>% 
  ts_plot()

fit1 <- Arima(test.ar_5$y, order = c(1,0,0), xreg = test.ar_5$X)
fit2 <- Arima(test.ar_20$y, order = c(1,0,0), xreg = test.ar_20$X)
fit3 <- Arima(test.ar_100$y, order = c(1,0,0), xreg = test.ar_100$X)
cor(fitted(fit1),test.ar_5$y )^2
cor(fitted(fit2),test.ar_20$y )^2
cor(fitted(fit3),test.ar_100$y )^2

# SETAR MODEL
test.setar_5 <- Model.setar(n_t[2], n_x[1], 
                            # err_sd = 0.6, 
                            alpha = 0.4, 
                            gamma = 0.8, 
                            th_r = -1, 
                            sigma = 2.5, 
                            delta = 3)

test.setar_20 <- Model.setar(n_t[2], n_x[2], 
                             # err_sd = 0.6, 
                             alpha = 0.4, 
                             gamma = 0.8, 
                             th_r = -1, 
                             sigma = 1, 
                             delta = 2)

test.setar_100 <- Model.setar(n_t[2], n_x[3], 
                              # err_sd = 0.6, 
                              alpha = 0.4, 
                              gamma = 0.8, 
                              th_r = -1, 
                              sigma = 1, 
                              delta = 2)

fit1 <- setar(test.setar_5$y, mL = 1, mH = 1, th = -1)
fit2 <- setar(test.setar_20$y, mL = 1, mH = 1, th = -1)
fit3 <- setar(test.setar_100$y, mL = 1, mH = 1, th = -1)
cor(fitted(fit1)[-1],test.setar_5$y[-1] )^2
cor(fitted(fit2)[-1],test.setar_20$y[-1] )^2
cor(fitted(fit3)[-1],test.setar_100$y[-1] )^2

ts_c(test.setar_5$y, 
     test.setar_20$y, 
     test.setar_100$y) %>% ts_plot

ts_c(`5 Xs` = (test.lin_5$y), 
     `5 Xs SetAR` = (test.setar_5$y)) %>% 
  ts_plot()

ts_c(`5 Xs AR` = (test.ar_5$y), 
     `5 Xs SetAR` = (test.setar_5$y)) %>% 
  ts_plot()

ts_c(`20 Xs` = (test.lin_20$y), 
     `20 Xs SetAR` = (test.setar_20$y)) %>% 
  ts_plot()

ts_c(`100 Xs` = test.lin_100$y, 
     `100 Xs SETAR` = test.setar_100$y) %>% 
  ts_plot()

# Markov-Switching Model

test.msw_5 <- Model.MS(n_t[2], n_x[1], 
                       rho = 0.9,
                       sigma_1 = 0.75,
                       sigma_2 = 0.5)

test.msw_20 <- Model.MS(n_t[2], n_x[2],  
                        rho = 0.9,
                        sigma_1 = 1.4,
                        sigma_2 = 0.5)

test.msw_100 <- Model.MS(n_t[2], n_x[3],  
                         rho = 0.9,
                         sigma_1 = 1,
                         sigma_2 = 0.5)

y <- test.msw_5$y
X <- test.msw_5$X
mod1 = lm(y~X)
mod.mswm1 <- msmFit(mod1,k=2,p=1,sw=rep(TRUE, n_x[1] + 3),control=list(parallel=TRUE))
summary(mod.mswm1)

y <- test.msw_20$y
X <- test.msw_20$X
mod2 = lm(y~X)
mod.mswm2 <- msmFit(mod2,k=2,p=1,sw=rep(TRUE, n_x[2] + 3),control=list(parallel=TRUE))
summary(mod.mswm2)

y <- test.msw_100$y
X <- test.msw_100$X
mod3 = lm(y~X)
summary(mod3)
mod.mswm3 <- msmFit(mod3,k=2,p=1,sw=rep(TRUE, n_x[3] + 3),control=list(parallel=TRUE))
summary(mod.mswm3)

ts_c(test.msw_5$y, 
     test.msw_20$y, 
     test.msw_100$y) %>% ts_plot

ts_c(`5 Xs` = (test.lin_5$y), 
     `5 Xs MSW` = (test.msw_5$y)) %>% 
  ts_plot()

ts_c(`5 Xs AR` = (test.ar_5$y), 
     `5 Xs MSWAR` = (test.msw_5$y)) %>% 
  ts_plot()

ts_c(`20 Xs` = (test.lin_20$y), 
     `20 Xs MSW` = (test.msw_20$y)) %>% 
  ts_plot()

ts_c(`100 Xs` = test.lin_100$y, 
     `100 Xs SETAR` = test.msw_100$y) %>% 
  ts_plot()

ts_c(`Linear` = test.lin_20$y,
     `AR` = test.ar_20$y,
     `SETAR` = test.setar_20$y,
     `MS` = test.msw_20$y) %>% 
  ts_plot()


# 5 Nonlinear Model

test.nonlin_5 <- Model.nonlinear(n_t[2], n_x[1], 
                       rho = 0.8,
                       sigma = 0.6)

test.nonlin_20 <- Model.nonlinear(n_t[2], n_x[2], 
                                  rho = 0.8,
                                  sigma = 0.6)

test.nonlin_100 <- Model.nonlinear(n_t[2], n_x[3], 
                                   rho = 0.8,
                                   sigma = 0.6)

ts_c(test.nonlin_5$y, 
     test.nonlin_20$y, 
     test.nonlin_100$y) %>% ts_plot


library(ggthemes)

fig1 <- ts_c(`Linear` = test.lin_20$y,
             `AR` = test.ar_20$y,
             `SETAR` = test.setar_20$y,
             `MS` = test.msw_20$y) %>%
  ts_tbl() %>% 
  ggplot() +
  geom_line(aes(y = value,
                x = time, 
                color = id),
  ) + 
  geom_hline(yintercept=0 , colour="#999999", size = .5, alpha = 0.5) +
  theme_tufte(base_size = 12, base_family = "serif", ticks = TRUE) +
    scale_color_brewer(palette="Dark2") + 
  labs(x = "Time",
       y = "Value")+
  theme(axis.text = element_text(size=10))+
  scale_y_continuous(
    expand = expand_scale(0.005),
    limits = c(-30, 40),
    breaks = seq(-30, 40, by = 10)) +
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  theme(legend.text=element_text(size=10),
        legend.position = "bottom",
        legend.title = element_blank())
fig1

main_dir <- "//adb.intra.admin.ch/Userhome$/SECO-01/U80839737/data/Documents/2_Forschung/4_Projects/_inflation_NC_HFdata/"

ggsave(paste0(main_dir,"/_results/","dgps.png"), 
       width = 8, height = 4.5, dpi = 300, units = "in", device = "png")


# Full sample regressions ----
library("texreg")

# Linear Model
lin_small <- Model.linear(n_t[1], n_x[1], err_sd = 1.1)
lin_large <- Model.linear(n_t[4], n_x[1], err_sd = 0.7)

fit_small <- lm(lin_small$y ~ lin_small$X)
names(fit_small$coefficients) <- c('intercept','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10')
cor(fitted(fit_small),lin_small$y )^2
fit_large <- lm(lin_large$y ~ lin_large$X)
names(fit_large$coefficients) <- c('intercept','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10')
cor(fitted(fit_large),lin_large$y )^2



mytable <- texreg(list(fit_small, fit_large), label = "tab:lin",
                  caption = "Linear Model",
                  float.pos = "h", return.string = TRUE, bold = 0.05, stars = 0,
                  custom.note = "Coefficients with $p < 0.05$ in \\textbf{bold}.",
                  digits = 3, leading.zero = FALSE, omit.coef = "Inter")

# ARIMAX Model
arima_small <- Model.arimax(n_t[1], n_x[1], rho = 0.8, err_sd_y = 3, err_sd_x = 1.15)
arima_large <- Model.arimax(n_t[4], n_x[1], rho = 0.8, err_sd_y = 3.8, err_sd_x = 1.35)

fit_small <- Arima(arima_small$y, order = c(1,0,0), xreg = arima_small$X)
fit_large <- Arima(arima_large$y, order = c(1,0,0), xreg = arima_large$X)
cor(fitted(fit_small), arima_small$y )^2
cor(fitted(fit_large), arima_large$y )^2

names(fit_small$coef) <- c('ar1','intercept','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10')
names(fit_large$coef) <- c('ar1','intercept','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10')

mytable <- texreg(list(fit_small, fit_large), label = "tab:lin",
                  caption = "AR(1)-X Model",
                  float.pos = "h", return.string = TRUE, bold = 0.05, stars = 0,
                  custom.note = "Coefficients with $p < 0.05$ in \\textbf{bold}.",
                  digits = 3, leading.zero = FALSE, omit.coef = "Inter")

# SETAR Model

setar_small <- Model.setar(n_t[1], n_x[1], 
                            # err_sd = 0.6, 
                            alpha = 0.4, 
                            gamma = 0.8, 
                            th_r = -1, 
                            sigma = 1.5, 
                            delta = 7.5)

setar_long <- Model.setar(n_t[4], n_x[2], 
                             # err_sd = 0.6, 
                             alpha = 0.4, 
                             gamma = 0.8, 
                             th_r = -1, 
                             sigma = 1.2, 
                             delta = 2.5)

fit_small <- setar(setar_small$y, mL = 1, mH = 1, th = -1)
fit_long <- setar(setar_long$y, mL = 1, mH = 1, th = -1)
cor(fitted(fit_small)[-1],setar_small$y[-1] )^2
cor(fitted(fit_long)[-1],setar_long$y[-1] )^2


