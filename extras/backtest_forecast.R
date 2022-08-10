library(RegTree)
library(data.table)

dt <- fread("/users/seth/data/price_fcasting_data_week_2022_08_05.csv")
names(dt)[1:2] <- c("ref_date", "tgt") 


dtt <- dt[weekdays(dt$ref_date) == "Tuesday"]
dtw <- dt[weekdays(dt$ref_date) == "Wednesday"]
dth <- dt[weekdays(dt$ref_date) == "Thursday"]

dates_t <- as.Date(dtt$ref_date)
dates_w <- as.Date(dtw$ref_date)
dates_h <- as.Date(dth$ref_date)

yt <- as.matrix(dtt[ , 2, with=FALSE])
Xt <- as.matrix(dtt[ , -c(1,2), with=FALSE])

yw <- as.matrix(dtw[ , 2, with=FALSE])
Xw <- as.matrix(dtw[ , -c(1,2), with=FALSE])

yh <- as.matrix(dth[ , 2, with=FALSE])
Xh <- as.matrix(dth[ , -c(1,2), with=FALSE])


out_t <- reg_forest(yt,Xt,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                  weight_by_mse = FALSE, weight_pow = 2, min_obs = 50)
out_w <- reg_forest(yw,Xw,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                    weight_by_mse = FALSE, weight_pow = 2, min_obs = 50)
out_h <- reg_forest(yh,Xh,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                    weight_by_mse = FALSE, weight_pow = 2, min_obs = 50)

sequoia <- c(out_t$Trees, out_t$Trees)
y_out <- StdFitField(Xw, sequoia)

Fit <- data.table("true_val" = c(yw), "fit" = c(y_out[[1]]))

pretty_plot(Fit)
mean(Fit$true_val^2, na.rm=TRUE)
mean((Fit$fit-Fit$true_val)^2, na.rm=TRUE)


mean(out$true_vals, na.rm=TRUE)
mean(out$out_of_sample)

pretty_plot(cbind(out$true_vals, out$out_of_sample))

mean(((out$out_of_sample) - y)^2, na.rm=TRUE)
mean(out$true_vals^2, na.rm=TRUE)

mean(sign(out$out_of_sample) == sign(out$true_vals), na.rm = TRUE)




# ---- Run a backtest over x days (this takes time) -----------------
tst_size <- 300
max_obs <- max(which(is.finite(dt$tgt)))
test_idx <- seq(max_obs-tst_size+1, max_obs)
test_dates <- dt$ref_date[test_idx]
y_fit <- rep(NA, tst_size) # j <- 1
y_true <- dt$tgt[test_idx]
for(j in seq(tst_size)){
  print(j)
  idx <- test_idx[j]
  dt_tst <- dt[seq(idx), ] # backtest dates
  dt_tst[(idx-4):idx, 2] <- NA 
  # one weekday at a time
  dtt <- dt_tst[weekdays(dt_tst$ref_date) == "Tuesday"]
  dtw <- dt_tst[weekdays(dt_tst$ref_date) == "Wednesday"]
  dth <- dt_tst[weekdays(dt_tst$ref_date) == "Thursday"]
  
  yt <- as.matrix(dtt[ , 2, with=FALSE])
  Xt <- as.matrix(dtt[ , -c(1,2), with=FALSE])
  
  yw <- as.matrix(dtw[ , 2, with=FALSE])
  Xw <- as.matrix(dtw[ , -c(1,2), with=FALSE])
  
  yh <- as.matrix(dth[ , 2, with=FALSE])
  Xh <- as.matrix(dth[ , -c(1,2), with=FALSE])
  
  out_t <- reg_forest(yt,Xt,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                      weight_by_mse = FALSE, weight_pow = 2, min_obs = 50)
  out_w <- reg_forest(yw,Xw,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                      weight_by_mse = FALSE, weight_pow = 2, min_obs = 50)
  out_h <- reg_forest(yh,Xh,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                      weight_by_mse = FALSE, weight_pow = 2, min_obs = 50)
  
  sequoia <- c(out_t$Trees, out_w$Trees, out_h$Trees)
  XX <- as.matrix(dt_tst[ , -c(1,2), with=FALSE])
  y_out <- StdFitField(as.matrix(dt_tst[ , -c(1,2), with=FALSE]), sequoia)
  y_fit[j] <- tail(y_out[[1]], 1)
}



# vs "out of bag"
# oob <- out$out_of_sample[test_idx]

# backtest_one_week <- cbind(y_fit, y_true, oob)
# write.csv(backtest_one_week, "/Users/seth/data/backtests/dg_week.csv")

pretty_plot(cbind(y_fit, y_true))

out$mean_abs_feature_contribution
colnames(X)

mean((y_fit-y_true)^2)
mean((y_true)^2)

sum(sign(y_fit) == sign(y_true))



# dt <- dt[seq(min(which(is.finite(dt$diff_pct_growth))), NROW(dt))]
# 
# names(dt)[1] <- "ref_date_zj40pq" # unique name
# dates <- as.Date(dt$ref_date_zj40pq)
# y <- dt$price_month_change
# X <- as.matrix(dt[ , c("pct_chng_ibes", "ibes_v_reported_yoy" ,
#                        "ibes_v_ibes_yoy", "pct_cum_growth",
#                        "pct_day_growth", "diff_pct_growth",
#                        "gap_v_consensus", "lag_price_pct_change",
#                        "lag_price_week_change",  "lag_price_month_change"), with=FALSE])


