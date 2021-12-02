
run_forecast <- function(tgt, dt, lib, regress = FALSE){
  MF <- process_MF(dt[series_name == tgt], dt[series_name != tgt], LHS_lags = 3, pub_date_name = NULL)
  X <- process(MF, lib, detrend = TRUE, pub_date_name = NULL)
  W <- dcast(X, ref_date ~ series_name, value.var = "value") # wide data
  tgt0 <- paste(tgt, "0")
  setcolorder(W, c("ref_date", tgt0))
  
  y <- unlist(W[ , tgt0, with = FALSE])
  dates <- W$ref_date
  setcolorder(W, c("ref_date", tgt0))
  W <- as.matrix(  W[ , -c(1,2), with = FALSE])
  
  # W[abs(W) > 4] <- NA # W is scaled, so this drops outliers
  if(regress){
    idx <- abs(W) > 5 | is.na(W)
    W[idx] <- 5*sign(W[idx]) # W is scaled, so this drops outliers
  } 
  
  out <- reg_forest(y, W, steps = 1, draws = 2000, regression = regress) # estimate model for one step ahead; pretty slow
  
  last_obs <- nrow(out$fit)
  raw_res <- data.table("ref_date" = dates[seq(last_obs)], "fit" = out$fit, "true" = out$true_vals)
  # pretty_plot(tail(raw_res, 48))
  
  df <- X[series_name == tgt0]
  
  # undo processing
  fit_1 <- df$standardize_scale[seq(last_obs)]*out$fit + df$standardize_center[seq(last_obs)] + df$low_frequency_trend[seq(last_obs)]
  fit_level <-exp(shift(df$level_value, 1)[seq(last_obs)]+fit_1)
  
  res <- data.table("series_name" = tgt, "ref_date" = dates[seq(last_obs)], "as_of" = Sys.Date(), "fit" = c(fit_1),  
                    "level_fit" = c(fit_level), "value" = Diff(df$level_value[seq(last_obs)]), "level_value" = df$level_value[seq(last_obs)])
  
  # pretty_plot(res[ , .(ref_date, fit, true_value)])
  
  out <- list("dt" = res,
              "raw" = out)
  return(out)
}

run_backtest <- function(tgt, dt, lib, regress = FALSE, periods = 24){
  as_of_dates <- c(as.Date(tail(dt[series_name == tgt]$pub_date, periods) - 2)) # include latest ?
  res <- data.table()
  # as_of_dates <- head(as_of_dates,6)
  
  for(j in seq(length(as_of_dates))){
    dte = as_of_dates[j]
    print(dte)
    MF <- process_MF(dt[series_name == tgt], dt[series_name != tgt], LHS_lags = 3, as_of = dte)
    X <- process(MF, lib, detrend = TRUE)
    W <- dcast(X, ref_date ~ series_name, value.var = "value") # wide data
    tgt0 <- paste(tgt, "0")
    setcolorder(W, c("ref_date", tgt0))
    
    y <- unlist(W[ , tgt0, with = FALSE])
    dates <- W$ref_date
    setcolorder(W, c("ref_date", tgt0))
    
    W <- as.matrix(  W[ , -c(1,2), with = FALSE])
    
    if(regress){
      idx <- abs(W) > 5 | is.na(W)
      W[idx] <- 5*sign(W[idx]) # W is scaled, so this drops outliers
    } 
    
    out <- reg_forest(y, W, steps = 1, draws = 2000, regression = regress, orthogonal = FALSE) # estimate model for one step ahead; pretty slow
    
    last_obs <- nrow(out$fit)
    raw_res <- data.table("ref_date" = dates[seq(last_obs)], "fit" = out$fit, "true" = out$true_vals)
    # pretty_plot(tail(raw_res, 24))
    
    df <- X[series_name == tgt0]
    
    # undo processing
    fit_1 <- df$standardize_scale[seq(last_obs)]*out$fit + df$standardize_center[seq(last_obs)] + df$low_frequency_trend[seq(last_obs)]
    fit_level <-exp(shift(df$level_value, 1)[seq(last_obs)]+fit_1)
    
    tmp <- data.table("ref_date" = dates[seq(last_obs)], "fit" = c(fit_1), "level_fit" = c(fit_level), "as_of" = dte)
    res <- rbind(res, tail(tmp,1))
  }
  
  tmp <- dt[series_name == tgt , .(ref_date, value)]
  setnames(tmp, "value", "level_value")
  tmp[ , value := Diff(log(value))]
  
  res <- merge(res, tmp, all.x = TRUE)
  res[ , series_name := tgt]
  setcolorder(res, c("series_name", "ref_date", "as_of", "fit", "level_fit", "value", "level_value"))
  # pretty_plot(res[ , .(ref_date, fit, pct_value)])
  
  
  return(res)
  
}
