# In sample estimates using a static data snapshot
library(data.table)
library(dateutils)
library(RegTree)
library(ranger)
library(randomForest)
library(macroeconomicdata)

# Load Data and library with instructions for adjustment
# load("./extras/final_lib.RData")
# load("./extras/final_dataset.RData")
# save.image("./extras/data_2022_05.RData")

load("./extras/data_2022_05.RData")

# only use SA import and export data
DT <- DT[!series_name%in%c("imports", "exports")]
LIB[series_name == "exports", series_name := "exports sa"]
LIB[series_name == "imports", series_name := "imports sa"]

#  Our in-sample results are from the day before CPI is published, i.e. they use
# the maximum plausible amount of data. The last publication date is
# 2022-01-13

# "consumer price index cpi" 
# "core consumer prices" 
# "personal consumption expenditures price index"
# "core pce price index"
# tgt <-  "core pce price index"

cpi_tgts <- c("consumer price index cpi", 
              "core consumer prices",
              "cpi com less food energy", 
              "cpi nondurables",
              "cpi durables", 
              "cpi food", 
              "cpi com", 
              "cpi serv", 
              "cpi serv less food energy")

pce_tgts <- c("core pce price index",
              "personal consumption expenditures price index")

best_rhs <- c("gasoline prices sa", 
              "gasoline prices weekly sa", 
              "ISM Manufacturing Prices", 
              "ISM Non Manufacturing Prices",
              "CRB Index",
              "T5YIE",
              "weekly economic index",
              "Brent", 
              "redbook index sa",
              "S&P GSCI",
              "inflation expectations 10y",
              "unemployment rate")

DT <- DT[series_name%in%c(cpi_tgts, pce_tgts, best_rhs)]
DT <- fill_pub_dates(DT) # belt and suspenders
DT <- DT[pub_date < as.Date("2022-01-12")]# 2022-01-12  # for CPI measures
# DT <- DT[ref_date >= as.Date("2010-01-01")]

# weekly economic index 0 
# 1.165390e-04                           9.595110e-05                           8.463063e-05 
# Brent 0                          Heating oil 0                     redbook index sa 0 
# 6.887176e-05                           6.731969e-05                           6.352745e-05 
# S&P GSCI 0           inflation expectations 10y 0                             Gasoline 0 
# 5.415007e-05                           3.510763e-05                           3.319780e-05 
# labor force participation rate 0                non manufacturing pmi 0                    manufacturing pmi 0 
# 3.311736e-05                           3.031027e-05                           2.899485e-05 
# Crude Oil 0                                 Gold 0

# OUT <- run_forecast(tgt = cpi_tgts[1],
#                     dt = DT,
#                     lib = LIB,
#                     detrend = FALSE,
#                     type = "standard",
#                     as_of_date = NULL)




out_df <- data.table("target"=cpi_tgts, "rf_full"=0, "rfrn_full"=0,
                     "rf_sub"=0, "rfrn_sub"=0)
stored_output <-list()

# DT[series_name == "core consumer prices"]


# DT <- DT[pub_date < as.Date("2021-12-23")] # for PCE measures

# we will be forecasting December, and all OOB estimates will use data
# that would have been available on 2022-01-13 to predict cpi
# max(DT[series_name == "consumer price index cpi"]$ref_date) 

idx <- 1

for(j in 1:length(cpi_tgts)){
  for(k in c(1,2)){
    if(k==1) regress <- TRUE
    else regress <- FALSE
    print(paste(cpi_tgts[j], regress))
    OUT <- run_forecast(tgt = cpi_tgts[j],
                        dt = DT,
                        lib = LIB,
                        detrend = FALSE,
                        regress = regress,
                        as_of_date = NULL)
    x <- Diff(log(DT[series_name == cpi_tgts[j]]$value))
    dates <- DT[series_name == cpi_tgts[j]]$ref_date
    mod <- lm(x ~ shift(x))
    x_hat <- cbind(matrix(1, length(x),1), shift(x))%*%mod$coefficients
    sq_er_ar <- (x-x_hat)^2
    sq_er_rf <- (c(x,NA)-OUT$dt$fit)^2
    to_date <- which(dates >= as.Date('2013-10-01'))
    subsamp <- which(dates >= as.Date('2013-10-01') & dates < as.Date('2020-01-01'))
    
    mse_ar <- sqrt(mean(sq_er_ar[to_date], na.rm=TRUE))
    mse_rf <- sqrt(mean(sq_er_rf[to_date], na.rm=TRUE))
    full_samp <- mse_rf/mse_ar 
    mse_ar <- sqrt(mean(sq_er_ar[subsamp], na.rm=TRUE))
    mse_rf <- sqrt(mean(sq_er_rf[subsamp], na.rm=TRUE))
    sub_samp <- mse_rf/mse_ar 
    if(k==1){
      out_df[target==cpi_tgts[j], rfrn_full:=full_samp]
      out_df[target==cpi_tgts[j], rfrn_sub:=sub_samp]
    }else{
      out_df[target==cpi_tgts[j], rf_full:=full_samp]
      out_df[target==cpi_tgts[j], rf_sub:=sub_samp]
    }
    stored_output[[idx]] <- OUT$dt
    idx <- idx + 1
  }
}

out_df

# save.image(file = "./extras/InSampleSmall.RData")

# ------- Plots and Feature Contributions ------------------------

tgt <- "cpi serv less food energy"

OUT <- run_forecast(tgt = tgt,
                    dt = DT,
                    lib = LIB,
                    detrend = FALSE,
                    regress = TRUE,
                    as_of_date = NULL)

# make a pretty plot

W <- dcast(OUT$X, ref_date ~ series_name, value.var = "value") # wide data
colnames(W)


par(family = "Times New Roman", mar = c(8.1,4.1,4.1,2.1) )
dt_plot <- OUT$dt[ref_date >= as.Date("2013-09-30")]

# png(filename ="./extras/cpi_insample_fit3.png", width = 600, height = 400)

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


# dev.off()

# ------------------- Feature Contributions --------------------------------

fc <- OUT$feature_contribution[ref_date > as.Date("2013-09-30")]

mfc <- colMeans(abs(fc[,-1, with=FALSE]), na.rm = TRUE)

sort(mfc, decreasing = TRUE)

LIB[series_name == "unemployed persons"]

ts.plot(DT[series_name == "weekly economic index"]$value)

fc_plot <- fc[ , .(ref_date, `gasoline prices sa 0`,  
                   `gasoline prices weekly sa 0`, 
                   `CRB Index 0`,
                   `T5YIE 0`,
                   `consumer price index cpi 1`)]

par(family = "Times New Roman", mar = c(5.1,4.1,4.1,10.2) )

# png(filename ="./extras/feature_contribution.png", width = 600, height = 400)

matplot(fc_plot$ref_date, fc_plot[,-1,with=FALSE], type = 'l', lty = c(1, 2,3,1,2), 
        col = c("black", "steelblue4", "steelblue3", "steelblue2", "steelblue"), lwd = 1, xlab = "Date", ylab = "Featrue Contribution", 
        bty = "n")
grid(col = "lightslategrey")
legend("right", c("Gasoline Monthly", "Gasoline Weekly", "CRB Index", "5 Year Inflation Expectations", "Lagged CPI"), 
       inset = c(-.33,0), xpd = TRUE,
       lty = c(1,2,3,1,2), lwd = 1,
       col = c("black", "steelblue4", "steelblue3", "steelblue2", "steelblue"), bty = "n",
       pt.cex=1, cex=.8)
title("Feature Contributions, Top 5 Variables", font.main=1)

# dev.off()

pretty_plot(fc[ , .(ref_date, `gasoline prices sa 0`,  
                    `gasoline prices weekly sa 0`, 
                    `CRB Index 0`,
                    `T5YIE 0`,
                    `ISM Non Manufacturing Prices 0`,
                    `ISM Manufacturing Prices 0`,
                    `Pickup Truck 0`)])




