library(RegTree)
library(data.table)
library(dateutils)

load("/Users/seth/data/inflation_nowcast_paper/final_dataset.RData")
load("/Users/seth/data/inflation_nowcast_paper/final_lib.RData")

DT[series_name == "consumer price index cpi"]

DT <- DT[series_name != "inflation rate mom"] # drop this one (publication_date is before cpi)

# feature contributions depend on our "as_of_date". I've set it here to be the 
# day before the last publication date, so we're using all available data.
# Obviously for an earlier date, if a series is not available, its feature 
# contribution will be zero :-(. 
out <- run_forecast(tgt = "consumer price index cpi", dt = DT, lib = LIB, 
                    as_of_date = "2022-05-10", type = "regression")  # type = "standard"

# just for fun
pretty_plot(out$dt[ , .(ref_date, fit, value)])

# which feature contributions matter most?
sort(out$raw$mean_abs_feature_contribution, decreasing = TRUE)[1:10]

# plot the key feature contributions
pretty_plot(tail(out$feature_contribution[ , .(ref_date, 
                                          `gasoline prices sa 0`, 
                                          `consumer price index cpi 1`, 
                                          `ism serv prices 0`, 
                                          `CRB Index 0`,
                                          `empire state manu fut price paid 0`,
                                          `bank lending rate 0`,
                                          `inflation expectations 10y 0`)],96))

fc <- out$feature_contribution[ , .(ref_date, 
                                    `gasoline prices sa 0`, 
                                    `consumer price index cpi 1`, 
                                    `ism serv prices 0`, 
                                    `CRB Index 0`,
                                    `empire state manu fut price paid 0`,
                                    `bank lending rate 0`,
                                    `inflation expectations 10y 0`)]

write.csv(fc, file = "/Users/seth/data/inflation_nowcast_paper/feature_contribution_regression.csv", row.names = FALSE)
