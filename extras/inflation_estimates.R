
rm(list=ls()) # clear workspace

library(data.table)
library(seasonal) # Christoph's package
library(dateutils) # My own package
library(RegTree)
library(randomForest)
library(ranger)
library(macroeconomicdata)
library(jsonlite)
library(httr)

windowsFonts(
  A=windowsFont("Arial Black"),
  B=windowsFont("Bookman Old Style"),
  C=windowsFont("Comic Sans MS"),
  D=windowsFont("Corbel")
)
par(family = "D")

DT <- fread("C:/Users/seton/Dropbox/macro/usa_data_sa.csv") # load data
lib <- fread("C:/Users/seton/Dropbox/macro/library.csv")
lib <- lib[country == "united states"]
DT[ , ref_date := as.Date(ref_date)]
DT[ , pub_date := as.Date(pub_date)]

DT <- DT[!series_name%in%c("non farm payrolls sa", "crude oil rigs sa")]

# --- for early plots in slides ----------------------------------------------------

fred_series <- c("CPIAUCSL", "PCEPI")
inflation <-lapply(fred_series, FUN = Get_FRED_Data)  
inflation <- rbindlist(inflation)
inflation <- dcast(inflation, ref_date ~ series_name)
names(inflation) <- c("ref_date", "CPI level", "PCE level")
inflation[ , CPI := (1+pct_chng(`CPI level`))^12-1]
inflation[ , PCE := (1+pct_chng(`PCE level`))^12-1]
pretty_plot(x = inflation$ref_date, X = as.matrix(inflation[, .(CPI, PCE)]), title = "Inflation Measures")

dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/pce_cpi.png", width = 1400, height = 900, res = 200)

vs <- tail(inflation[ , .(ref_date, `CPI level`, CPI)], 12)
vs[ , `CPI level` := 3*`CPI level`/vs$`CPI level`[1]]
tmp <- vs$`CPI level`
ig <- 1.04^(1/12)
ig <- ig^seq(6)
tmp[7:12] <- tmp[6]*ig

X <- cbind(tmp - 3, vs$`CPI level`-3, vs$CPI)

matplot(vs$ref_date, X, type = 'l', lty = c(2,1,1), col = c("black", "black", "deeppink"), lwd = 2, xlab = "Date")
grid(col = "lightslategrey")
legend("topleft", c("4% Inflation in Levels", "CPI level", "CPI percent change"),
       col = c("black", "black", "deeppink"), lty = c(2,1,1), lwd = 2) #  bty = "n"
title("Inflation vs. Price Level")

dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/cpi_inf_v_lev.png", width = 1400, height = 900, res = 175)

keep_small <- c("consumer price index cpi",
                "Brent", "gasoline prices weekly",
                "corporate 10y t bill spread")
dt <- DT[series_name%in%keep_small]
dt <- dcast(dt, ref_date ~ series_name, value.var = "value")
dt <- dt[ref_date >= as.Date("2021-08-31")]
write.csv(dt, file = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/MF_example.csv", row.names = FALSE, na = " ")

# ----------------------------------------------------------------------------

# load commodity prices
CM <- fread("C:/Users/seton/Dropbox/market/commod_data.csv") # load data
CM[ , ref_date := as.Date(ref_date)]
clib <- data.table(series_name = unique(CM$series_name), needs_SA = FALSE, 
                   take_logs = TRUE, take_diffs = TRUE)

# library for all data
LIB <- rbind(lib[ , .(series_name, needs_SA, take_logs, take_diffs)], clib)
LIB[ , needs_SA := FALSE]

median(DT[series_name == "consumer price index cpi"]$pub_lag) # find pub_lag
selected <- unique(DT[pub_lag<=13]$series_name) # CPI is lagged 13 days from reference date

DT <- DT[series_name%in%selected & ref_date >= as.Date("1992-01-01")] # begin in 1992
CM <- CM[ref_date >= as.Date("1992-01-01")]

# rbind data sets
CM <- CM[ , .(country, series_name, ref_date, close)]
setnames(CM, "close", "value")
CM[ , pub_date := ref_date+1]
CM[ , pub_lag := 1]

DT <- rbind(DT, CM) # one data set
DT[ , country := NULL] # drop country since we're only dealing with the US

# --- add EA data -------------------------
ea_data <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\blended.csv")
ea_data[ , ref_date := as.Date(ref_date)]
ea_data[ , pub_date := ref_date]
ea_data[ , pub_lag := 0]

DT <- rbind(DT, ea_data)

# ------ add EA library ----------------

ea_lib <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\hack_lib.csv")
LIB <- rbind(LIB, ea_lib)

# ----------- Backtest the model ------------------------
back_dates <- seq.Date(from = as.Date("2019-01-06"), to = Sys.Date(), by = "week") # 2021 to date
Fit <- matrix(NA, length(back_dates)+1, 2) # one and two step ahead forecasts
fcast_date <- rep(as.Date(NA), length(back_dates)) # ref_date for one step ahead forecast
for(j in seq(length(back_dates))){
  as_of <- back_dates[j] # current as of date
  print(as_of)
  dat <- process_MF(LHS = DT[series_name == "consumer price index cpi"], RHS = DT[series_name != "consumer price index cpi"],
                    LHS_lags = 3, RHS_lags = 0, as_of = as_of, frq = "month") # aggregate
  dt <- process(dat, lib = LIB) # make stationary
  cpi <- dt[series_name == "consumer price index cpi 0"] # LHS variable
  y <- cpi$value
  X <- dcast(dt[series_name != "consumer price index cpi 0"], ref_date ~ series_name, value.var = "value")
  dates <- X$ref_date
  X <- as.matrix(X[,-1,drop=FALSE]) # RHS variables
  
  out <- reg_forest(y, X, steps = 1) # estimate model for one step ahead
  # undo processing
  last_obs <- nrow(out$fit)
  fit_1 <- cpi$standardize_scale[seq(last_obs)]*out$fit + cpi$standardize_center[seq(last_obs)] + cpi$low_frequency_trend[seq(last_obs)]
  fcast_date[j] <- dates[last_obs] # date of one step ahead forecast
  
  out <- reg_forest(y, X, steps = 2) # estimate model for two steps ahead
  # undo processing
  last_obs <- nrow(out$fit)
  fit_2 <- cpi$standardize_scale[seq(last_obs)]*out$fit + cpi$standardize_center[seq(last_obs)] + cpi$low_frequency_trend[seq(last_obs)]
  
  Fit[j,1] <- tail(fit_1,1)
  Fit[j+1,2] <- tail(fit_2,1)
}

fcast_dates <- c(fcast_date, end_of_period(tail(fcast_date,1), period = "month", shift = 1))
Fit <- data.table("ref_date" = fcast_dates, Fit)
true_val <- DT[series_name == "consumer price index cpi"] # get true values
true_val[ , value := Diff(log(value))]
res <- merge(Fit, true_val[ , .(ref_date, value)], by = "ref_date", all.x = TRUE)

names(res) <- c("ref_date", "one step ahead", "two steps ahead", "true values")
res <- data.table("as_of" = c(back_dates, tail(back_dates,1)+7), res)

Res_old <- fread("~/inflation-nowcasting/init_results.csv")
res <- merge(res, Res_old[ , .(as_of, cleveland_fed)], by = "as_of", all = TRUE)
res[ , `true values` := 100*`true values`]
res[ , `one step ahead` := 100*`one step ahead`]
res[ , `two steps ahead` := 100*`two steps ahead`]
setcolorder(res, c("as_of", "ref_date", "true values", "one step ahead", "two steps ahead", "cleveland_fed"))

pretty_plot(X = as.matrix(res[,-c(1,2, 6),with=FALSE]), x = res$as_of, title = "Out of Sample Backtest, CPI")
dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/cpi_backtest.png", width = 1400, height = 900, res = 175)

pretty_plot(X = as.matrix(res[,-c(1,2,5),with=FALSE]), x = res$as_of, title = "Out of Sample Backtest, CPI")
dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/cpi_backtest_fed.png", width = 1400, height = 900, res = 175)

mean((res$`one step ahead` -  res$`true values`)^2, na.rm = TRUE)
mean((res$`cleveland_fed` -  res$`true values`)^2, na.rm = TRUE)

# Save results
write.csv(res, file = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/backtest_results.csv", row.names = FALSE, na = " ")

# -------------------------------------------------------

# --- Prediction with current dataset ----------------------

dat <- process_MF(LHS = DT[series_name == "consumer price index cpi"], RHS = DT[series_name != "consumer price index cpi"],
                  LHS_lags = 3, RHS_lags = 0, as_of = NULL, frq = "month") # aggregate
dt <- process(dat, lib = LIB) # make stationary
cpi <- dt[series_name == "consumer price index cpi 0"]
y <- cpi$value
X <- dcast(dt[series_name != "consumer price index cpi 0"], ref_date ~ series_name, value.var = "value")
dates <- X$ref_date
X <- as.matrix(X[,-1,drop=FALSE])

out <- reg_forest(y, X, steps = 1) # estimate model
# undo processing
last_obs <- nrow(out$fit)
fit <- cpi$standardize_scale[seq(last_obs)]*out$fit + cpi$standardize_center[seq(last_obs)] + cpi$low_frequency_trend[seq(last_obs)]

# plot against true values
fit <- data.table("ref_date" = dates[seq(last_obs)], "fitted cpi" = c(fit))
true_val <- DT[series_name == "consumer price index cpi"]
true_val[ , value := Diff(log(value))]
res <- merge(true_val[ , .(ref_date, value)], fit, by = "ref_date", all = TRUE)
setnames(res, "value", "true value")
res[ , `true value`:= 100*`true value`]
res[ , `fitted cpi` := 100*`fitted cpi`]
pretty_plot(tail(res, 48), ylab = "Monthly Percent Change", title = "CPI Nowcast")

dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/cpi_september.png", width = 1400, height = 900, res = 175)


tail(res)
tmp <- DT[series_name == "consumer price index cpi"]
(1.00378388)*tail(tmp$value,1)

# Analysis: which series contribute the most?


first_split(out$Trees, X_names = colnames(X))

all_splits(out$Trees, X_names = colnames(X))

# ------------- Short sample analysis -------------------------------------

DT <- DT[ref_date >= as.Date("2010-01-01")]

dat <- process_MF(LHS = DT[series_name == "consumer price index cpi"], RHS = DT[series_name != "consumer price index cpi"],
                  LHS_lags = 3, RHS_lags = 0, as_of = NULL, frq = "month") # aggregate
dt <- process(dat, lib = LIB) # make stationary
cpi <- dt[series_name == "consumer price index cpi 0"]
y <- cpi$value
X <- dcast(dt[series_name != "consumer price index cpi 0"], ref_date ~ series_name, value.var = "value")
dates <- X$ref_date
X <- as.matrix(X[,-1,drop=FALSE])

out <- reg_forest(y, X, steps = 1) # estimate model
# undo processing
last_obs <- nrow(out$fit)
fit <- cpi$standardize_scale[seq(last_obs)]*out$fit + cpi$standardize_center[seq(last_obs)] + cpi$low_frequency_trend[seq(last_obs)]

# plot against true values
fit <- data.table("ref_date" = dates[seq(last_obs)], "fitted cpi" = c(fit))
true_val <- DT[series_name == "consumer price index cpi"]
true_val[ , value := Diff(log(value))]
res <- merge(true_val[ , .(ref_date, value)], fit, by = "ref_date", all = TRUE)
setnames(res, "value", "true value")
res[ , `true value`:= 100*`true value`]
res[ , `fitted cpi` := 100*`fitted cpi`]
pretty_plot(tail(res, 48), ylab = "Monthly Percent Change", title = "CPI Nowcast")

# pretty_plot(tail(res, 48), ylab = " ", title = " ", xlab = " ")

dev.print(png, filename = "C:/Users/seton/Dropbox/Quantagon/inflation_nowcasting/cpi_september_short.png", width = 1400, height = 900, res = 175)

first_split(out$Trees, X_names = colnames(X))

all_splits(out$Trees, X_names = colnames(X))



















res <- cbind(out$true_vals, out$fit)
colnames(res) <- c("true vals", "fitted vals")
pretty_plot(res)




first_obs <- function(x) min(which(is.finite(x)))
fst_obs <- sapply(X, FUN = first_obs)






