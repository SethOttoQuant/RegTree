
rm(list=ls()) # clear workspace

library(data.table)
library(seasonal) # Christoph's package
library(dateutils) # My own package, on CRAN
library(RegTree)
library(randomForest)
library(ranger)
library(macroeconomicdata) # My own packages
library(jsonlite)
library(httr)
library(tradingeconomics)
login("068185ACDE82499:864AE58966F24E9")

# DT <- fread("C:/Users/seton/Dropbox/macro/usa_data_sa.csv") # load data
# lib <- fread("C:/Users/seton/Dropbox/macro/library.csv")

DT <- fread("~/Dropbox/macro/usa_data_sa.csv") # load data
lib <- fread("~/Dropbox/macro/library.csv")

lib <- lib[country == "united states"]
lib <- lib[ , .(series_name, frequency, needs_SA, take_logs, take_diffs)]

DT[ , ref_date := as.Date(ref_date)]
DT[ , pub_date := as.Date(pub_date)]
DT[ , country := NULL]

DT <- DT[!series_name%in%c("non farm payrolls sa", "crude oil rigs sa")]

# ------- Add data not already in DT --------------------

TEnew <- c("MBA Mortgage Market Index", "MBA Mortgage Refinance Index")
frednew <- c("TOTBKCR", "RELACBW027SBOG", "CLSACBW027SBOG", "FRGSHPUSM649NCIS", "FRGEXPUSM649NCIS")
# all comercial debt, real estate loans, consmer loans, Cass freight shipments (NSA), Cass freigh expenditures (NSA) 

# Additional data from Trading Economics
new_data <- lapply(TEnew, FUN = TE_historical_data, country = "United States", initDate = "1990-01-01") 
dat <- rbindlist(new_data)
dat[ , pub_lag := 5]
dat[ , pub_date := ref_date + pub_lag]
dat[ , country := NULL]
setcolorder(dat, names(DT))
# Additional data from FRED
frd <- lapply(frednew, FUN = Get_FRED_Data, observation_start = "1990-01-01")
frd <- rbindlist(frd)
frd[series_name%in%c("TOTBKCR", "RELACBW027SBOG", "CLSACBW027SBOG") , pub_lag := 5]
frd[series_name%in%c("FRGSHPUSM649NCIS", "FRGEXPUSM649NCIS"), pub_lag := 12]
frd[ , pub_date := ref_date + pub_lag]
frd[series_name%in%c("FRGSHPUSM649NCIS", "FRGEXPUSM649NCIS"), value := log(value)]
# the two Cass series need seasonal adjustment
frd[series_name == "FRGSHPUSM649NCIS", value := sa_inplace(value, ref_date)]
frd[series_name == "FRGEXPUSM649NCIS", value := sa_inplace(value, ref_date)]
setcolorder(frd, names(DT))

dt <- rbind(DT, dat)
dt <- rbind(dt,frd)

libdat <- auto_library(dat)
# libdat[ , frequency:= NULL]
# tmp <- dcast(dat, ref_date ~ series_name, value.var = "value")
# pretty_plot(tmp, legend_pos = "topleft")

libfrd <- auto_library(frd) # Fine
# libfrd[ , frequency:= NULL]

lib <- rbind(lib, libdat)
lib <- rbind(lib, libfrd)

# ----------------------------------------------------------------------------

# load commodity prices
CM <- fread("~/Dropbox/market/commod_data.csv") # load data
CM[ , ref_date := as.Date(ref_date)]
clib <- data.table(series_name = unique(CM$series_name), frequency = "daily",
                   needs_SA = FALSE,  take_logs = TRUE, take_diffs = TRUE)

# library for all data
LIB <- rbind(lib, clib)
LIB[ , needs_SA := FALSE]

median(DT[series_name == "consumer price index cpi"]$pub_lag) # find pub_lag
selected <- unique(DT[pub_lag<=13]$series_name) # CPI is lagged 13 days from reference date

DT <- DT[series_name%in%selected & ref_date >= as.Date("1992-01-01")] # begin in 1992
CM <- CM[ref_date >= as.Date("1992-01-01")]

# rbind data sets
CM <- CM[ , .(series_name, ref_date, close)]
setnames(CM, "close", "value")
CM[ , pub_date := ref_date+1]
CM[ , pub_lag := 1]

DT <- rbind(DT, CM) # one data set

DT <- DT[!(series_name%in%c("gasoline prices sa", "gasoline prices", "gasoline prices weekly sa" ))]
unique(DT$series_name)

# --- Prediction with current dataset ----------------------

dat <- process_MF(LHS = DT[series_name == "consumer price index cpi"], RHS = DT[series_name != "consumer price index cpi"],
                  LHS_lags = 3, RHS_lags = 1, as_of = NULL, frq = "month") # aggregate as.Date("2021-11-09")
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
gas <- dt[series_name == "gasoline prices weekly 0"]
gas[ , `gas price` := value/2]
res <- merge(res, gas[, .(ref_date, `gas price`)], all.x = TRUE)
pretty_plot(tail(res, 48), ylab = "Monthly Percent Change", title = "CPI Nowcast")

# dev.print(png, filename = "~/Dropbox/system2/inflation/cpi_november.png", width = 1400, height = 900, res = 150)

tail(res)


(1+tail(res$`fitted cpi`,1)/100)^12
lev <- cpi$level_value
exp(lev[max(which(is.finite(lev)))] + tail(res$`fitted cpi`,1))

# Analysis: which series contribute the most?

first_split(out$Trees, X_names = colnames(X))

all_splits(out$Trees, X_names = colnames(X))


