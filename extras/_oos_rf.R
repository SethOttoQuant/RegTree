rm(list=ls()) # clear workspace

library(data.table)
library(jsonlite)
library(httr)
library(tidyverse)
library(readxl)

library(randomForest)
library(ranger)

library(tsbox)
library(seasonal) # Christoph's package

library(dateutils) # My own package
library(RegTree)
library(macroeconomicdata)
library(database)

setwd(here::here())
setwd("../")
main_dir<- getwd()

refresh_data <- FALSE

run_forecast <- function(tgt_id, 
                         tgt_lags = 3, 
                         dt, 
                         lib, 
                         detrend = FALSE, 
                         regress = FALSE, 
                         as_of_date = NULL, 
                         hsteps = 1){
  
  # Define mixed-frequency data matirx
  dat.MF <- process_MF(LHS = dt[series_name == tgt_id],
                       RHS = dt[series_name != tgt_id],
                       LHS_lags = tgt_lags,
                       RHS_lags = 0,
                       as_of = as_of_date,
                       frq = "month") # aggregate
  
  # Make data stationary
  X <- process(dat.MF,
               lib = lib,
               detrend = detrend,
               pub_date_name = NULL)   
  
  W <- dcast(X, 
             ref_date ~ series_name, 
             value.var = "value") # wide data
  
  tgt0 <- paste(tgt_id, "0")
  setcolorder(W, c("ref_date", tgt0)) # this should be the case anyhow
  
  y <- W[[tgt0]] 
  dates <- W$ref_date
  setcolorder(W, c("ref_date", tgt0))
  W <- as.matrix(  W[ , -c(1,2), with = FALSE])
  
  if(regress){
    idx <- abs(W) > 5 | is.na(W)
    W[idx] <- 5*sign(W[idx]) # W is scaled, so this drops outliers
  } 
  
  # estimate model for h-step ahead; pretty slow
  out <- reg_forest(y, W, 
                    steps = hsteps, 
                    draws = 2000, 
                    regression = regress) 
  
  last_obs <- nrow(out$fit)
  raw_res <- data.table("ref_date" = dates[seq(last_obs)], 
                        "fit" = out$fit, 
                        "true" = out$true_vals)
  
  df <- merge(data.table("ref_date"=dates), 
              X[series_name == tgt0], 
              by = "ref_date", 
              all.x = TRUE)
  
  df[ , value := Diff(level_value)] # not standardized/scaled
  fit <- out$out_of_sample
  upper_bound <- fit + sqrt(out$mse)
  lower_bound <- fit - sqrt(out$mse)
  
  # undo processing
  if(detrend){
    fit <- df$standardize_scale[seq(last_obs)]*fit + df$standardize_center[seq(last_obs)] + df$low_frequency_trend[seq(last_obs)]
    upper_bound <- df$standardize_scale[seq(last_obs)]*upper_bound + df$standardize_center[seq(last_obs)] + df$low_frequency_trend[seq(last_obs)]
    lower_bound <- df$standardize_scale[seq(last_obs)]*lower_bound + df$standardize_center[seq(last_obs)] + df$low_frequency_trend[seq(last_obs)]
  }else{
    fit <- df$standardize_scale[seq(last_obs)]*fit + df$standardize_center[seq(last_obs)] 
    upper_bound <- df$standardize_scale[seq(last_obs)]*upper_bound + df$standardize_center[seq(last_obs)]
    lower_bound <- df$standardize_scale[seq(last_obs)]*lower_bound + df$standardize_center[seq(last_obs)]
  }
  
  if(lib[series_name == tgt_id]$take_diffs){
    fit_level <-shift(df$level_value, 1)[seq(last_obs)]+fit
    upper_bound_level <-shift(df$level_value, 1)[seq(last_obs)]+upper_bound
    lower_bound_level <-shift(df$level_value, 1)[seq(last_obs)]+lower_bound
  }else{
    fit_level <- fit
    upper_bound_level <- upper_bound
    lower_bound_level <- lower_bound
  }
  if(lib[series_name == tgt_id]$take_logs){
    fit_level <- exp(fit_level)
    upper_bound_level <- exp(upper_bound_level)
    lower_bound_level <- exp(lower_bound_level)
    df[ , level_value := exp(level_value)]
  }
  
  res <- data.table("series_name" = tgt_id, 
                    "ref_date" = dates[seq(last_obs)], 
                    "as_of" = as_of_date, 
                    "fit" = c(fit), 
                    "upper_bound" = upper_bound, 
                    "lower_bound" = lower_bound, 
                    "level_fit" = c(fit_level), 
                    "upper_bound_level" = upper_bound_level,
                    "lower_bound_level" = lower_bound_level,
                    "value" = df$value[seq(last_obs)], 
                    "level_value" = df$level_value[seq(last_obs)])
  
  out <- list("dt" = res,
              "raw" = out,
              "X" = X)
  return(out)
}


if (isTRUE(refresh_data)){
  #---- 1) Establish the data set ------------------------------------------------
  
  # a) General US Macro data -----------------------------------------------------
  MD <- fread(paste0(main_dir, "/_data/usa_data_sa.csv")) 
  MD[ , ref_date := as.Date(ref_date)]
  MD[ , pub_date := as.Date(pub_date)]
  
  # drop series --> doubles and nsa where sa is available
  MD <- MD[!series_name%in%c("non farm payrolls sa", 
                             "crude oil rigs sa", 
                             "inflation expectations", 
                             "gasoline prices weekly", 
                             "average weekly hours",
                             "challenger job cuts",
                             "fiscal expenditure",
                             "gasoline prices",
                             "government revenues",
                             "natural gas stocks change",
                             "redbook index")]
  MD[ , country := NULL] # drop country since we're only dealing with the US
  
  lib <- fread(paste0(main_dir, "/_data/library.csv"))
  lib <- lib[country == "united states"]
  
  # b) Commodity prices ----------------------------------------------------------
  CM <- fread(paste0(main_dir, "/_data/commod_data.csv"))
  CM[ , ref_date := as.Date(ref_date)]
  CM <- CM[ , .(country, series_name, ref_date, close)]
  setnames(CM, "close", "value")
  CM[ , pub_date := ref_date+1]
  CM[ , pub_lag := 1]
  CM[ , country := NULL] # drop country since we're only dealing with the US
  
  clib <- data.table(series_name = unique(CM$series_name), 
                     needs_SA = FALSE, 
                     take_logs = TRUE, 
                     take_diffs = TRUE)
  
  # c) Additional FRED data ------------------------------------------------------
  Fred <- rbindlist(lapply(c("T10YIE", "T5YIE"), 
                           FUN = Get_FRED_Data, 
                           observation_start = "1992-01-01"))
  Fred[ , pub_date := ref_date]
  Fred[ , pub_lag := 0]
  
  FredLib <- data.table("series_name" = c("T10YIE", "T5YIE"), 
                        "needs_SA" = 0, 
                        "take_logs" = 0, 
                        "take_diffs" = 0 )
  
  # d) Additional Macrobond data -------------------------------------------------
  mb_ids_d <- c("uslead0497",	
                "ovx",	
                "ctotusd",	
                "cesiusd",
                "wocaes0029",	
                "usepundnewsindex",
                "usepundeqmuncindex",	
                "uslead0019",	
                "ussurv02742")
  
  mb_ids_w <- c("ussurv6195",	
                "usprod0685",	
                "ustran0032",	
                "usrigcountbasinw0064",	
                "ussurv0387",
                "ussurv01117",	
                "ussurv03506",	
                "fiber0001",	
                "usbust3493")
  
  meta <- read_xlsx(paste0(main_dir, "/_data/macrobond.xlsx"), sheet = "lib") %>% 
    rename(`id` = "mb_id")
  
  mb_d <- macrobond(mb_ids_d, class = "tbl")
  mb_w <- macrobond(mb_ids_w, class = "tbl") 
  
  dta_mb <- bind_rows(mb_d, mb_w) %>% 
    left_join(meta, by = "id") %>% 
    rename(`ref_date` = "time") %>% 
    mutate(pub_date = ref_date + pub_lag)
  
  MB <- dta_mb[, .(series_name, ref_date, value, pub_date, pub_lag)]
  
  MB_lib <- dta_mb[, .(series_name, needs_SA,take_logs, take_diffs)] %>% 
    distinct(series_name, .keep_all = TRUE)
  
  # e) Merge data sets -----------------------------------------------------------
  DT <- rbind(MD, 
              CM, 
              MB, 
              Fred) # one data set
  
  # f) Establish library for all data --------------------------------------------
  LIB <- rbind(lib[ , .(series_name, needs_SA, take_logs, take_diffs)], 
               clib, 
               FredLib, 
               MB_lib)
  LIB[ , needs_SA := FALSE]
  
  
  # g) A first overview: ---------------------------------------------------------
  DT %>% 
    select(id = series_name, 
           time = ref_date, 
           value) %>% 
    ts_tbl() %>% 
    ts_summary()
  
  DT <- DT[ref_date >= as.Date("1992-01-01")] # begin in 1992 --> why is that 
  
  save(DT, file = paste0(main_dir, "/_results/final_dataset.Rdata"))
  save(LIB, file = paste0(main_dir, "/_results/final_lib.Rdata"))
}
#---- 2) Prepare pseudo real-time dataset --------------------------------------
load(paste0(main_dir, "/_results/final_dataset.Rdata"))
load(paste0(main_dir, "/_results/final_lib.Rdata"))

# ----------- Backtest the model ------------------------

# CPI -------------------------------------------------------------------------- 

tgt_id <- "consumer price index cpi"
tgt_lags <- 3

pub_lag_target <- median(DT[series_name == tgt_id]$pub_lag) # find pub_lag
selected <- unique(DT[pub_lag <= pub_lag_target]$series_name) # CPI is lagged 13 days from reference date
DATA <- DT[series_name%in%selected]


# Part 1
back_dates <-  seq(as.Date("2019-01-05"),  #as.Date("2013-10-05"), 
                   to = as.Date("2021-12-31"), 
                   by = "week")

Results <-list()

for(j in seq(length(back_dates[1:145]))){
  as_of <- back_dates[j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.cpi_part1.Rdata"))

Results <-list()

for(j in seq(length(back_dates[146:290]))){
  as_of <- back_dates[145 + j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.cpi_part2.Rdata"))

Results <-list()

for(j in seq(length(back_dates[291:430]))){
  as_of <- back_dates[290+j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.cpi_part3.Rdata"))

# Core CPI ---------------------------------------------------------------------

tgt_id <- "core consumer prices"
tgt_lags <- 3

pub_lag_target <- median(DT[series_name == tgt_id]$pub_lag) # find pub_lag
selected <- unique(DT[pub_lag <= pub_lag_target]$series_name) # CPI is lagged 13 days from reference date
DATA <- DT[series_name%in%selected]

# Part 1
back_dates <-  seq(as.Date("2013-10-05"),  
                   to = as.Date("2021-12-31"), 
                   by = "week")

Results <-list()

for(j in seq(length(back_dates[1:145]))){
  as_of <- back_dates[j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.cpi_core_part1.Rdata"))

Results <-list()

for(j in seq(length(back_dates[146:290]))){
  as_of <- back_dates[145 + j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.cpi_core_part2.Rdata"))

Results <-list()

for(j in seq(length(back_dates[291:430]))){
  as_of <- back_dates[290+j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.cpi_core_part3.Rdata"))

#  PCE -------------------------------------------------------------------------

tgt_id <- "personal consumption expenditures price index"
tgt_lags <- 3

pub_lag_target <- median(DT[series_name == tgt_id]$pub_lag) # find pub_lag
selected <- unique(DT[pub_lag <= pub_lag_target]$series_name) # CPI is lagged 13 days from reference date
DATA <- DT[series_name%in%selected]

# Part 1
back_dates <-  seq(as.Date("2013-10-05"),  
                   to = as.Date("2021-12-31"), 
                   by = "week")

Results <-list()

for(j in seq(length(back_dates[1:145]))){
  as_of <- back_dates[j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.pce_part1.Rdata"))

Results <-list()

for(j in seq(length(back_dates[146:290]))){
  as_of <- back_dates[145 + j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.pce_part2.Rdata"))

Results <-list()

for(j in seq(length(back_dates[291:430]))){
  as_of <- back_dates[290+j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.pce_part3.Rdata"))


#  PCE CORE --------------------------------------------------------------------

tgt_id <- "core pce price index"
tgt_lags <- 3

pub_lag_target <- median(DT[series_name == tgt_id]$pub_lag) # find pub_lag
selected <- unique(DT[pub_lag <= pub_lag_target]$series_name) # CPI is lagged 13 days from reference date
DATA <- DT[series_name%in%selected]

# Part 1
back_dates <-  seq(as.Date("2013-10-05"),  
                   to = as.Date("2021-12-31"), 
                   by = "week")

Results <-list()

for(j in seq(length(back_dates[1:145]))){
  as_of <- back_dates[j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.pce_core_part1.Rdata"))

Results <-list()

for(j in seq(length(back_dates[146:290]))){
  as_of <- back_dates[145 + j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.pce_core_part2.Rdata"))

Results <-list()

for(j in seq(length(back_dates[291:430]))){
  as_of <- back_dates[290+j] # current as of date
  print(as_of)
  
  # RANDOM FOREST
  out_rf <- run_forecast(tgt_id, 
                         tgt_lags = 3,
                         dt = DATA, 
                         lib = LIB, 
                         detrend = FALSE,
                         regress = FALSE,
                         as_of_date = as_of, 
                         hsteps = 1)
  
  # RANDOM FOREST RANDOM NODE
  out_rfrn <- run_forecast(tgt_id, 
                           tgt_lags = 3,
                           dt = DATA, 
                           lib = LIB, 
                           detrend = FALSE,
                           regress = TRUE,
                           as_of_date = as_of, 
                           hsteps = 1)
  
  out <- list(rf = out_rf, 
              rfrn = out_rfrn)
  
  Results[[length(Results)+1]] = out
}

save(Results, file = paste0(main_dir, "/_results/Results.pce_core_part3.Rdata"))


# Expand sample size
# Define the four different targets
# Expand forecasting horizon
# Define the two different models

