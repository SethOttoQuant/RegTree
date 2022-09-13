library(data.table)
library(dateutils)
library(RegTree)


library(RegTree)
library(data.table)

dt <- fread("/users/seth/data/price_fcasting_data_2022_08_24.csv")
names(dt)[1:2] <- c("ref_date", "tgt") # date and target varible

# In sample over-fitting example
dates <- dt$ref_date
y <- dt$tgt
X <- as.matrix(dt[,-c(1,2), with=FALSE])
out_is <- reg_forest(y,X,draws=1000, type = "standard", return_trees = FALSE, steps = NULL, 
                    weight_by_mse = TRUE, weight_pow = 2) # in sample example

dt_is <- data.table("ref_date" = dates[out_is$idx],  "true_value" = c(out_is$true_vals),
                    "in_sample_fit" = c(out_is$in_samp_fit))

pretty_plot(tail(dt_is, 300))
pretty_plot(dt_is)


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

# x <- matrix(rnorm(500), 100, 5)
# y <- x%*%c(1,3,2,1,-1) + rnorm(100)
# bestsplit(x, y, seq(1,100), 3, 10)
# out <- fast_cut(x[,1], y, 10)
# out[4]

out_t <- reg_forest(yt,Xt,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                    weight_by_mse = FALSE, weight_pow = 2, max_obs = 50, min_obs = 10)
out_w <- reg_forest(yw,Xw,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                    weight_by_mse = FALSE, weight_pow = 2, max_obs = 50, min_obs = 10)
out_h <- reg_forest(yh,Xh,draws=1000, type = "standard", return_trees = TRUE, steps = NULL, 
                    weight_by_mse = FALSE, weight_pow = 2, max_obs = 50, min_obs = 10)


dt_t <- data.table("true_value" = c(out_t$true_vals),
                    "fit" = c(out_t$in_samp_fit))

pretty_plot(dt_t)

mean(sign(dt_t$true_value) == sign(dt_t$fit), na.rm = TRUE)
1 - mean((dt_t$true_value - dt_t$fit)^2, na.rm = TRUE)/mean(dt_t$true_value^2, na.rm=TRUE)

bob <- data.table(out_t$Trees[[10]])

bob[ , variance := V7/V8]

yy <- scale(yt)
ts.plot(yy)

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