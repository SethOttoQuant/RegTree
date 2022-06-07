# In sample estimates using a static data snapshot
library(data.table)
library(dateutils)
library(RegTree)
library(ranger)
library(randomForest)
library(macroeconomicdata)
library(tradingeconomics)
login("068185ACDE82499:864AE58966F24E9")

macro_path <- "/tmp/usa_data_sa.csv" 
commod_path <- "/tmp/commod_data.csv"
library_path <- "/tmp/library.csv"
cars_path <- "/tmp/cars_latest.csv"

# get data
MD <- fread(macro_path) # can also get this via Dropbox API
MD[ , ref_date := as.Date(ref_date)]
MD[ , pub_date := as.Date(pub_date)]
MD <- MD[ref_date >= as.Date("1992-01-01")]
MD[ , country := NULL]
MD[ , pub_date := NULL]
MD[ , pub_lag := NULL]

CM <- fread(commod_path) # load data: saved locally
CM <- CM[ , .(series_name, ref_date, close)]
setnames(CM, "close", "value")
CM[ , ref_date := as.Date(ref_date)]
CM <- CM[ref_date >= as.Date("1992-01-01")]

# lib <- fread("C:/Users/seton/Dropbox/macro/library.csv")
lib <- fread(library_path)
lib <- lib[country == "united states"]
lib <- lib[ , .(series_name, needs_SA, take_logs, take_diffs)]

clib <- data.table(series_name = unique(CM$series_name), 
                   needs_SA = FALSE, 
                   take_logs = TRUE, 
                   take_diffs = TRUE)

FredDay <- rbindlist(lapply(c("T10YIE", "T5YIE"), 
                            FUN = Get_FRED_Data, 
                            observation_start = "1992-01-01"))
FredDay <- FredDay[ , .(series_name, ref_date, value)]
FredDay[ , ref_date := as.Date(ref_date)]

# CPIFABSL  - food and beverage
# CPITRNSL  - transport
# CPIHOSSL  - housing
# CUSR0000SETA02  - used cars and trucks
# CPIMEDSL - medical
# CPIEDUSL - education
# CPIRECSL  - recreation
# CPIAPPSL  - apparel
# CPIOGSSL  - other


fred_inflation <- c("CPIOGSSL", "CPIEDUSL", "CPIMEDSL", "CPITRNSL", "CPIHOSSL", "CPIFABSL", "CPIRECSL", "CPIAPPSL")
fred_monthly <- c(fred_inflation, "CUSR0000SETA02")

FredMonth <- rbindlist(lapply(fred_monthly, 
                              FUN = Get_FRED_Data, 
                              observation_start = "1992-01-01")
)


FredMonth <- FredMonth[ , .(series_name, ref_date, value)]
FredMonth[ , ref_date := end_of_period(as.Date(ref_date))]

FredLib <- data.table("series_name" = c("T10YIE", "T5YIE", fred_monthly), 
                      "needs_SA" = FALSE, 
                      "take_logs" = FALSE, 
                      "take_diffs" = FALSE)
FredLib[series_name%in%fred_monthly, take_logs:=TRUE]
FredLib[series_name%in%fred_monthly, take_diffs:=TRUE]

# # Look at the actual series we are dealing with
# yy <- DT[series_name == "personal consumption expenditures durable goods" ]
# pretty_plot(yy[ , .(ref_date, value)])

# mpl <- DT[ , median(pub_lag), by = series_name]
# mpl <- mpl[V1 <= 30] # pub lags < 30 (i.e. median pub lag of pce durable goods)
# DT <- DT[series_name%in%mpl$series_name] # keep only stuff that is timely

# Additional data from Trading Economics
TEnew <- c("ISM Manufacturing Prices", "ISM Non Manufacturing Prices") # "MBA Mortgage Refinance Index"
new_data <- lapply(TEnew, FUN = TE_historical_data, country = "United States", initDate = "1992-01-01") 
TE <- rbindlist(new_data)
TE[ , country := NULL]
setcolorder(TE, names(MD))
TElib <- auto_library(TE) # fine
TElib[ , frequency:=NULL]

# These data to not seem very useful... 
# instock <- fread(instock_path)
# instock[ , series_name := paste(series_name, "instock", sep="_")]
# newin <- fread(newin_path)
# newin[ , series_name := paste(series_name, "newin", sep="_")]
# AC <- rbind(instock, newin)
# AClib <- data.table("series_name" = unique(AC$series_name), 
#                     needs_SA = FALSE, take_logs = FALSE, take_diffs = TRUE)
# AC <- setcolorder(AC, c("series_name", "ref_date"))
# AC[ , ref_date := as.Date(ref_date)]

CG <- fread(cars_path)
setcolorder(CG, c("label", "value_date", "value"))
names(CG) <- c("series_name", "ref_date", "value")
CG[ , ref_date := as.Date(ref_date)]
CGlib <- data.table("series_name" = unique(CG$series_name), "needs_SA" = FALSE, 
                    "take_logs" = TRUE, "take_diffs" = TRUE)
# pretty_plot(DT[series_name=="gasoline prices weekly sa", .(ref_date, value)])

TE[ , pub_date := ref_date + 5]
TE[ , pub_lag := 5]
write.csv(TE, file = "/tmp/te.csv", row.names = FALSE)

FredDay[ , pub_date := ref_date + 1]
FredDay[ , pub_lag := 1]
write.csv(FredDay, file = "/tmp/fred_day.csv", row.names = FALSE)

write.csv(TElib, file = "/tmp/te_day.csv", row.names = FALSE)
write.csv(FredLib[1:2], file = "/tmp/fred_lib.csv", row.names = FALSE)


DT <- rbind(MD, CM, TE, FredDay, FredMonth, CG)
LIB <- rbind(lib, clib, FredLib, TElib, CGlib)
LIB[, needs_SA:=FALSE]

DT <- DT[!series_name%in%c("non farm payrolls sa", 
                           "crude oil rigs sa", 
                           "inflation expectations", 
                           "gasoline prices weekly", 
                           "average weekly hours",
                           "challenger job cuts",
                           "fiscal expenditure",
                           "gasoline prices",
                           "government revenues",
                           "natural gas stocks change",
                           "gasoline stocks change",
                           "redbook index")]

# consumer price index cpi

tgt <- "consumer price index cpi"

OUT <- run_forecast(tgt = tgt,
                    dt = DT,
                    lib = LIB,
                    detrend = FALSE,
                    regress = TRUE,
                    as_of_date = NULL)

# make a pretty plot

par(family = "Times New Roman", mar = c(8.1,4.1,4.1,2.1) )
dt_plot <- OUT$dt[ref_date >= as.Date("2013-09-30")]

matplot(dt_plot$ref_date, dt_plot[,.(value, fit)], type = 'l', lty = c(1, 2), 
        col = c("steelblue4" ,"black"), lwd = 1, xlab = "Date", ylab = "CPI", 
        bty = "n")
polygon(x = c(dt_plot$ref_date, rev(dt_plot$ref_date)),
        y = c(dt_plot$lower_bound, 
              rev(dt_plot$upper_bound)),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
grid(col = "lightslategrey")
legend("bottom", c("True CPI", "MO-RFRN"), inset = c(0,-.4), xpd = TRUE,
       lty = c(1,2), lwd = 1, col = c("steelblue4" ,"black"), bty = "n", horiz = TRUE)
title("Out of Bag Fit, All Available Data", font.main=1)

mse = mean((dt_plot$value - dt_plot$fit)^2)
mst = mean(dt_plot)

# save plot to disk

dev.print(png, "./extras/cpi_insample_fit2.png", width = 600, height = 400)

# ----------- comp with AR(1) model -----------------------------
dates <- DT[series_name == tgt]$dates
x <- Diff(log(DT[series_name == tgt]$value))
mod <- lm(x ~ shift(x))
x_hat <-cbind(matrix(1,length(x),1),shift(x))%*%mod$coefficients # don't mess up indexes!!
mse_ar <- sqrt(mean((x-x_hat)^2, na.rm = TRUE))
# ar <- data.table("ref_date" = dates, "ar_fit" = x_hat)
mse_rf <- 

# min(OUT$dt$ref_date), min(OUT$dt$value, na.rm = TRUE)

# ------------------- Feature Contributions --------------------------------

fc <- OUT$feature_contribution[ref_date > as.Date("2013-09-30")]

mfc <- colMeans(abs(fc[,-1, with=FALSE]), na.rm = TRUE)

sorted_mfc <- sort(mfc, decreasing = TRUE)

top_5 <- names(sorted_mfc)[1:5]
fc_5 <- fc[colnames(fc)%in%top_5]

OUT$raw$mean_abs_feature_contribution

