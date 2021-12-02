# Load needed libraries
library(RegTree)
library(dateutils)
library(data.table)
library(macroeconomicdata)
library(tradingeconomics)
login("068185ACDE82499:864AE58966F24E9")

# get data
DT <- fread("/tmp/usa_data_sa.csv") # can also get this via Dropbox API
# DT <- fread("C:/Users/seton/Dropbox/macro/usa_data_sa.csv") # pc version
DT[ , ref_date := as.Date(ref_date)]
DT[ , pub_date := as.Date(pub_date)]
DT <- DT[ref_date >= as.Date("1990-01-01")]
DT[ , country := NULL]
DT[ , pub_date := NULL]
DT[ , pub_lag := NULL]

CM <- fread("/tmp/commod_data.csv") # load data: saved locally

# lib <- fread("C:/Users/seton/Dropbox/macro/library.csv")
lib <- fread("/tmp/library.csv")
lib <- lib[country == "united states"]
lib <- lib[ , .(series_name, frequency, needs_SA, take_logs, take_diffs)]

# # Look at the actual series we are dealing with
# yy <- DT[series_name == "personal consumption expenditures durable goods" ]
# pretty_plot(yy[ , .(ref_date, value)])

# mpl <- DT[ , median(pub_lag), by = series_name]
# mpl <- mpl[V1 <= 30] # pub lags < 30 (i.e. median pub lag of pce durable goods)
# DT <- DT[series_name%in%mpl$series_name] # keep only stuff that is timely

TEnew <- c("MBA Mortgage Market Index") # "MBA Mortgage Refinance Index"
frednew <- c("TOTBKCR", "RELACBW027SBOG", "CLSACBW027SBOG", "FRGSHPUSM649NCIS", 
             "FRGEXPUSM649NCIS", "DFDHRC1Q027SBEA", "DMOTRC1Q027SBEA", "DREQRC1Q027SBEA",
             "PCEND", "DODGRC1Q027SBEA", "RETAILIRSA", "R423IRM163SCEN")
# all comercial debt, real estate loans, consumer loans, Cass freight shipments (NSA), Cass freigh expenditures (NSA) 

# Additional data from Trading Economics
new_data <- lapply(TEnew, FUN = TE_historical_data, country = "United States", initDate = "1990-01-01") 
dat <- rbindlist(new_data)
dat[ , country := NULL]
setcolorder(dat, names(DT))
# Additional data from FRED
frd <- lapply(frednew, FUN = Get_FRED_Data, observation_start = "1990-01-01")
frd <- rbindlist(frd)
frd[series_name%in%c("FRGSHPUSM649NCIS", "FRGEXPUSM649NCIS", "PCEND", "RETAILIRSA", "R423IRM163SCEN"), ref_date := end_of_period(ref_date)] # monthly
frd[series_name%in%c("DFDHRC1Q027SBEA", "DMOTRC1Q027SBEA", "DREQRC1Q027SBEA", "DODGRC1Q027SBEA"), ref_date := end_of_period(ref_date, period = "quarter")]
# Seasonally adjust
frd[series_name == "FRGSHPUSM649NCIS", value := sa_inplace(value, ref_date)]
frd[series_name == "FRGEXPUSM649NCIS", value := sa_inplace(value, ref_date)]
# reset names:
frd[series_name == "DFDHRC1Q027SBEA", series_name := "pce furniture and home goods"]
frd[series_name == "DMOTRC1Q027SBEA", series_name := "pce motor vehicles"]
frd[series_name == "DREQRC1Q027SBEA", series_name := "pce recreation"]
frd[series_name == "DODGRC1Q027SBEA", series_name := "pce durable other"]
frd[series_name == "PCEND", series_name := "pce non durable"]
frd[series_name=="RETAILIRSA", series_name := "retail inventory sales ratio"]
frd[series_name=="R423IRM163SCEN", series_name := "wholesale durable inventory sales ratio"]

setcolorder(frd, names(DT))

dt <- rbind(DT, dat)
dt <- rbind(dt,frd)

libdat <- auto_library(dat) # fine
# tmp <- dcast(dat, ref_date ~ series_name, value.var = "value")
# pretty_plot(tmp, legend_pos = "topleft")

libfrd <- auto_library(frd) # Fine
clib <- data.table(series_name = unique(CM$series_name), frequency = "daily",
                   needs_SA = FALSE,  take_logs = TRUE, take_diffs = TRUE) # lib for commodity prices

lib <- rbind(lib, libdat)
lib <- rbind(lib, libfrd)
lib <- rbind(lib, clib)

# Read in Yodlee data: both aggregated and disaggregated

rev <- fread("/tmp/daily.csv")
rev[ , total_amount := NULL]
names(rev) <- c("ref_date", "series_name", "value")
rev[ , ref_date := as.Date(ref_date)]
rev <- rev[order(series_name, ref_date)]
setcolorder(rev, c("series_name", "ref_date", "value"))


# agrev <- rev[ , sum(value), by = ref_date]
# setnames(agrev, "V1", "value")
# agrev[ , series_name := "yodlee"]
# agrev[ , pub_lag := 10]
# agrev[ , pub_date := ref_date + pub_lag]
# setcolorder(agrev, c("series_name", "ref_date", "value", "pub_date", "pub_lag"))

# library file

lib[, needs_SA:=FALSE]
revlib <- auto_library(rev)
revlib[ , take_logs := TRUE]
revlib[ , take_diffs := TRUE]
revlib[ , needs_SA := TRUE]
# revlib <- rbind(revlib, data.table("series_name" = "yodlee", "frequency" = "week", 
#                          "needs_SA" = TRUE, "take_logs" = TRUE, "take_diffs" = TRUE))
lib <- rbind(lib, revlib)

dt <- dt[series_name != "exports sa"]
dt <- dt[series_name != "imports sa"]

# dt <- rbind(dt, rev)

# save.image("~/tmp/data_backup.RData")
# load("~/tmp/data_backup.RData")

# pretty_plot(dt[series_name == "economic optimism index"]$value, dt[series_name == "economic optimism index"]$ref_date)



selected_pced <- c("personal consumption expenditures durable goods", "personal consumption expenditures services", "advance retail sales",
              "MBA Mortgage Market Index", "CLSACBW027SBOG", "retail inventory sales ratio", "pce furniture and home goods", "pce motor vehicles", "pce recreation",
              "pce durable other", "pce non durable", "composite pmi", "consumer confidence",
              "corporate 10y t bill spread", "mortgage applications", "redbook index sa", "non farm payrolls",
              "services pmi", "retail inventory sales ratio", "wholesale durable inventory sales ratio",
              "car production", "car registrations sa", "durable goods orders ex defense sa", "stock market")

selected_quarterly <- c(selected_pced, "wholesale inventories", "building permits",
                        "factory orders ex transportation", "wholesale durable inventory sales ratio" )

# ----- Nowcast durable good consumption -------------------------

pced <- run_forecast(tgt="personal consumption expenditures durable goods", rbind(dt[series_name%in%selected_pced], rev), lib, regress = TRUE)

# ------- Furniture and home goods ---------------------------------

pcefhg <- run_forecast(tgt="pce furniture and home goods", rbind(dt[series_name%in%selected_quarterly], rev), lib, regress = TRUE)
# pretty_plot(pcefhg$dt[ , .(ref_date, fit, value)])

# ------ PCE Motor Vehicles -------------------------------------------

pcemv <- run_forecast(tgt="pce motor vehicles", rbind(dt[series_name%in%selected_quarterly], rev), lib, regress = TRUE)
# pretty_plot(pcemv$dt[ , .(ref_date, fit, value)])

out <- rbind(tail(pced$dt, 1), tail(pcefhg$dt,1), tail(pcemv$dt,1))
write.csv(out, "/tmp/pred.csv", row.names = FALSE)

# And other fun stuff it's probably worth keeping track of:

dtout <- dt[series_name%in%c("pce recreation", "pce non durable", "pce durable other",
                             "retail inventory sales ratio", "wholesale durable inventory sales ratio",
                             "pce recreation", "personal consumption expenditures services",
                             "personal consumption expenditures durable goods", "pce furniture and home goods",
                             "pce motor vehicles")]

dtout[ , pct_change := pct_chng(value)]
dtout[ , max_date := max(ref_date), by = series_name]
dtout <- dtout[ref_date == max_date]
dtout[, max_date := NULL]

write.csv(out, "/tmp/latest.csv", row.names = FALSE)



