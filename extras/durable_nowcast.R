# Load needed libraries
library(RegTree)
library(dateutils)
library(data.table)
library(macroeconomicdata)
library(tradingeconomics)
login("068185ACDE82499:864AE58966F24E9")

args <- commandArgs(trailingOnly = TRUE)

print(args)

yodlee_path <- args[[1]]
macro_path <- args[[2]]
commod_path <- args[[3]]
library_path <- args[[4]]

# yodlee_path = "/tmp/daily.csv" # path to yodlee data
# macro_path = "~/dropbox/macro/usa_data_sa.csv" # path to US macro data
# commod_path = "~/dropbox/market/commod_data.csv" # path to commodity data
# library_path = "~/dropbox/macro/library.csv"

# get data
DT <- fread(macro_path) # can also get this via Dropbox API
# jim <- fread("~/Dropbox/market/commod_data.csv") # pc version
DT[ , ref_date := as.Date(ref_date)]
DT[ , pub_date := as.Date(pub_date)]
DT <- DT[ref_date >= as.Date("1990-01-01")]
DT[ , country := NULL]
DT[ , pub_date := NULL]
DT[ , pub_lag := NULL]

# CM <- fread(commod_path) # load data: saved locally

# lib <- fread("C:/Users/seton/Dropbox/macro/library.csv")
lib <- fread(library_path)
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
             "PCEND", "DODGRC1Q027SBEA", "RETAILIRSA", "R423IRM163SCEN",
             "R4232IM163SCEN")
# all comercial debt, real estate loans, consumer loans, Cass freight shipments (NSA), Cass freigh expenditures (NSA) 

# Additional data from Trading Economics
# new_data <- lapply(TEnew, FUN = TE_historical_data, country = "United States", initDate = "1990-01-01") 
# dat <- rbindlist(new_data)
# dat[ , country := NULL]
# setcolorder(dat, names(DT))

# Additional data from FRED
frd <- lapply(frednew, FUN = Get_FRED_Data, observation_start = "1990-01-01")
frd <- rbindlist(frd)
frd[series_name%in%c("FRGSHPUSM649NCIS", "FRGEXPUSM649NCIS", "PCEND", "RETAILIRSA", "R423IRM163SCEN", "R4232IM163SCEN"), ref_date := end_of_period(ref_date)] # monthly
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
frd[series_name=="R4232IM163SCEN", series_name := "wholesale furniture home goods inventory sales ratio"]

setcolorder(frd, names(DT))

dt <- rbind(DT,frd)
# dt <- rbind(dt, dat)
# libdat <- auto_library(dat) # fine
libfrd <- auto_library(frd) # Fine
# clib <- data.table(series_name = unique(CM$series_name), frequency = "daily",
#                    needs_SA = FALSE,  take_logs = TRUE, take_diffs = TRUE) # lib for commodity prices

# lib <- rbind(lib, libdat)
lib <- rbind(lib, libfrd)
# lib <- rbind(lib, clib)

# Read in Yodlee data: both aggregated and disaggregated
print("Loading Yodlee data")
rev <- fread(yodlee_path)
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
                        "factory orders ex transportation", "wholesale durable inventory sales ratio",
                        "wholesale furniture home goods inventory sales ratio")

# ----- Nowcast durable good consumption -------------------------
print("Estimating PCE Durable")
pced <- run_forecast(tgt="personal consumption expenditures durable goods", rbind(dt[series_name%in%selected_pced], rev), lib, regress = FALSE)

# ------- Furniture and home goods ---------------------------------
print("Estimating PCE furniture and home goods")
pcefhg <- run_forecast(tgt="pce furniture and home goods", rbind(dt[series_name%in%selected_quarterly], rev), lib, regress = FALSE)
# pretty_plot(pcefhg$dt[ , .(ref_date, fit, value)])

# ------ PCE Motor Vehicles -------------------------------------------
print("Estimating PCE motor vehicles")
pcemv <- run_forecast(tgt="pce motor vehicles", rbind(dt[series_name%in%selected_quarterly], rev), lib, regress = FALSE)
# pretty_plot(pcemv$dt[ , .(ref_date, fit, value)])

print("Writting output")
out <- rbind(tail(pced$dt, 1), tail(pcefhg$dt,1), tail(pcemv$dt,1))
write.csv(out[ , .(series_name, ref_date, as_of, fit, level_fit, value, level_value)], "/tmp/pred.csv", row.names = FALSE)

# And other fun stuff it's probably worth keeping track of:

print("Adding other cool stuff")
dtout <- dt[series_name%in%c("pce recreation", "pce non durable", "pce durable other",
                             "retail inventory sales ratio", "wholesale durable inventory sales ratio",
                             "pce recreation", "personal consumption expenditures services",
                             "personal consumption expenditures durable goods", "pce furniture and home goods",
                             "pce motor vehicles", "wholesale furniture home goods inventory sales ratio")]

dtout[ , level_value := value]
dtout[ , value := Diff(log(value))]
dtout[ , max_date := max(ref_date), by = series_name]
dtout <- dtout[ref_date == max_date] # just output the latest data
dtout[, max_date := NULL]


write.csv(dtout, "/tmp/latest.csv", row.names = FALSE)







