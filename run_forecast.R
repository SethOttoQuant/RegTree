library(RegTree)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

print(args)

file_path <- args[[1]]
if(length(args)>1){
  type <- args[[2]]
}else{
  type <- "standard"
}
if(length(args)>2){
  weight_pow <- args[[3]] # number. Default is 2
  if(weight_pow==0){  # set to zero to not weight
    weight_by_mse <- FALSE
  }else{
    weight_by_mse <- TRUE
  }
}else{
  weight_by_mse <- TRUE
  weight_pow<-2
}
if(length(args)>3){
  draws <- args[[3]]
}else{
  draws <- 1000
}

file_path <- "/users/seth/data/price_fcasting_data_2022_08_02.csv"
dt <- fread(file_path)

names(dt)[1] <- "ref_date_zj40pq"
dates <- as.Date(dt$ref_date_zj40pq)
y <- as.matrix(dt[ , 2, with=FALSE])
#  pct_chng_mean_eps, pct_chng_mean_revenue, lag_price_pct_change, lag_price_week_change,
X <- as.matrix(dt[ , -c(1,2), with=FALSE])

out <- reg_forest(y,X,draws=draws, type = type, return_trees = FALSE, steps = NULL, 
                  weight_by_mse = weight_by_mse, weight_pow = weight_pow)

Fit <- data.table("ref_date" = dates[out$idx], "true_vals" = c(out$true_vals), "fitted_vals" = c(out$out_of_sample), "in_sample_fit" = c(out$in_samp_fit))

write.csv(Fit, file = file_path, row.names = FALSE)






