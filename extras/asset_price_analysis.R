library(data.table)
library(dateutils)
library(RegTree)

dt <- fread("/users/seth/data/price_fcasting_data_2022_08_02.csv")
dt <- dt[publication_date>=min(dt[is.finite(diff_pct_growth)]$publication_date)]
names(dt)[1] <- "ref_date"

# pretty_plot(dt[ , .(ref_date, diff_pct_growth, price_pct_change)])

dates <- as.Date(dt$price_date)
y <- dt$price_twoweek_change
#  pct_chng_mean_eps, pct_chng_mean_revenue, lag_price_pct_change, lag_price_week_change,
X <- as.matrix(dt[ , .(pct_chng_ibes, ibes_v_reported_yoy, ibes_v_ibes_yoy, 
                       pct_day_growth, diff_pct_growth, gap_v_consensus, 
                       lag_price_pct_change, lag_price_week_change, lag_price_month_change)])
# y <- dt$diff_pct_growth
# n <- length(y)-5
# y <- y[1:n]
# #  pct_chng_mean_eps, pct_chng_mean_revenue, lag_price_pct_change, lag_price_week_change,
# X <- as.matrix(dt[1:n , .(intraday, price_pct_change, price_twoday_change, price_threeday_change, price_week_change)])


# out <- reg_forest(y,X,draws=2000, type = "alt", geom_par = .5, min_obs = 100)

out <- reg_forest(y, X, weight_by_mse = TRUE)

Fit <- data.table("true_vals" = c(out$true_vals), "fitted_vals" = c(out$out_of_sample))
Fit <- tail(Fit, 400)
r2 <- 1 - sum((out$true_vals -out$out_of_sample)^2, na.rm=TRUE) / sum((out$true_vals - mean(out$true_vals, na.rm = TRUE))^2, na.rm = TRUE)
r2

pretty_plot(Fit)

insamp <- data.table("true_vals" = c(out$true_vals), "in_sample" = c(out$in_samp_fit))
pretty_plot(tail(insamp, 200))

var(y, na.rm = TRUE)
out$mse

out$mean_abs_feature_contribution
colnames(X)

oos <- out$out_of_sample




names(dt)