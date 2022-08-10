library(data.table)
library(dateutils)
library(RegTree)


library(RegTree)
library(data.table)

dt <- fread("/users/seth/data/price_fcasting_data_week_2022_08_05.csv")
names(dt)[1] <- "ref_date" # unique name
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