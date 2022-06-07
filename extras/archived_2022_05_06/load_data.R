library(data.table)
library(seasonal) # Christoph's package
library(dateutils) # My own package
library(RegTree)
library(randomForest)
library(ranger)
library(macroeconomicdata)
library(jsonlite)
library(httr)

# ------ Causality Link -----------------------------------
CL <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\CausalityLink\\causality20210928.csv")
CL[ , V1 := NULL]
CL_small <- CL[conceptkey == "usa-inflation-aggregated_score_top10"]
CL_small[ , link_weight := NULL]
CL_small[ , rel_link_weight := NULL]
CL_small[ , sourcekey_rnk := NULL]
CL_small[ , conceptkey := NULL]
CL_small[ , date := as.Date(date)]
setnames(CL_small, "date", "ref_date")

# ------- LinkUp -----------------------------------------
LU <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\LinkUp\\ea_inflation_sets.csv")
LU[ , ref_date := as.Date(date, format = "%m/%d/%y")]
LU <- LU[ref_date >= as.Date("2007-12-01")]
LU[ , idx := us_n_created - us_n_deleted]
LU_small <- LU[ , .(ref_date, idx)]

LU <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\LinkUp\\ea_duration_update.csv")
LU[ , ref_date := as.Date(day, format = "%m/%d/%y")]
LU_small <- merge(LU_small, LU[ , .(ref_date, unique_active_job_count)], all = TRUE)

# ------- Revelio ---------------------------------------
construction <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\Revelio\\construction_posting_summary.csv")
construction[ , post_week := as.Date(post_week)]
names(construction) <- c("ref_date", "build_avg_time_to_fill", "build_avg_salary", "build_n_postings")

# food <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\Revelio\\eaglealpha_food_long_file.csv")
# food[ , post_week := as.Date(post_week)]
# names(food) <- c("date", "food_avg_time_to_fill", "food_avg_salary", "food_n_postings")

restaurant <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\Revelio\\restaurant_posting_summary.csv")
restaurant[ , post_week := as.Date(post_week)]
names(restaurant) <- c("ref_date", "restaurant_avg_time_to_fill", "restaurant_avg_salary", "restaurant_n_postings")

RL <- merge(construction, restaurant, by = "ref_date", all = TRUE)

# ------- SpaceKnow ---------------------------------------

car1 <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\SpaceKnow\\light_vehicle_manufacturing_asi_usa.csv")
car1[ , Date := as.Date(Date)]
names(car1) <- c("ref_date", "car_asi")

car2 <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\SpaceKnow\\light_vehicle_manufacturing_sai_usa.csv")
car2[ , Date := as.Date(substr(Date, 1, 8), format = "%m/%d/%y")]
names(car2) <- c("ref_date", "30_car_sai", "24_car_sai", "12_car_sai")

air <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\SpaceKnow\\cargo_airport_asi_usa.csv")
air[ , Date := as.Date(Date)]
names(air) <- c("ref_date", "air_asi")

SK <- merge(car1, car2, by = "ref_date", all = TRUE)
SK <- merge(SK, air, by = "ref_date", all = TRUE)

# ------- Ascential --------------------------

hack <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\blended.csv")
Asc <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\Ascential\\prices_20141229_20211003.csv")
Asc[ , V7 := NULL]
setnames(Asc, "START_DATE", "ref_date")
Asc[ , ref_date := as.Date(ref_date)]

obs <- dcast(Asc, ref_date ~  CATEGORY, value.var =  "NEWIN_PRODUCT_COUNT")
obs$Jumpsuits


sm <- rowSums(obs[,-1,with = FALSE], na.rm = TRUE)
mean(sm)

price <- dcast(Asc, ref_date ~ CATEGORY, value.var = "AVG_PRODUCT_PRICE")
finite <- all_finite(t(price[,-1,with=F]))
price <- price[ , c(TRUE, finite), with = FALSE]

newprice <- dcast(Asc, ref_date ~ CATEGORY, value.var = "AVG_PRODUCT_NEWIN_PRICE")
newprice <- newprice[ , c(TRUE, finite), with = FALSE]
newprice <- as.matrix(price[,-1,with=FALSE]) - as.matrix(newprice[,-1,with=FALSE])
newprice <- data.table("ref_date" = price$ref_date, newprice)

pidx <- data.table("ref_date" = price$ref_date, 
                   "price" = rowMeans(scale(as.matrix(price[,-1,with=FALSE])), na.rm = TRUE),
                   "newprice" = rowMeans(scale(as.matrix(newprice[,-1,with=FALSE])), na.rm = TRUE))

names(newprice)[-1] <- paste("new", names(newprice)[-1])

price_long <- melt(price, id.vars = "ref_date", variable.name = "series_name", na.rm = TRUE)

price_lib <- auto_library(price_long)
price_lib[ , frequency := NULL]
price_lib[ , needs_SA := TRUE]

newprice_long <- melt(newprice, id.vars = "ref_date", variable.name = "series_name", na.rm = TRUE)

newprice_lib <- auto_library(newprice_long)
newprice_lib[ , frequency := NULL]

pidx_long <- melt(pidx, id.vars = "ref_date", variable.name = "series_name", na.rm = TRUE)
pidx_lib <- auto_library(pidx_long)
pidx_lib[ , frequency := NULL]
pidx_lib[ , needs_SA := TRUE]


hack_lib <- fread("C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\hack_lib.csv")
hack_lib[ , take_logs := as.logical(take_logs)]
hack_lib[ , take_diffs := as.logical(take_diffs)]

hack[ , ref_date := as.Date(ref_date)]
Hack <- rbind(hack, pidx_long)

write.csv(Hack, file = "C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\blended_2.csv", row.names = FALSE, na = " ")

HACK <- rbind(Hack, price_long, newprice_long)
write.csv(HACK, file = "C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\blended_large.csv", row.names = FALSE, na = " ")

LIB <- rbind(hack_lib, price_lib, newprice_lib, pidx_lib)
write.csv(LIB, "C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\hack_lib_big.csv", row.names = FALSE)

tst <- newprice_long[series_name == "Jumpsuits"]

pretty_plot(X = as.matrix(tst$value), x = tst$ref_date)

# --- Merge it all together --------------
Hack <- merge(CL_small, LU_small, by = "ref_date", all = TRUE)
Hack <- merge(Hack, RL, by = "ref_date", all = TRUE)
Hack <- merge(Hack, SK, by = "ref_date", all = TRUE)

hack <- melt(Hack, id.vars = "ref_date", variable.name = "series_name", na.rm = TRUE)
write.csv(hack, file = "C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\blended.csv", row.names = FALSE, na = " ")


# ------- Creating a library file for the data -------

j <- 2

names(Hack)[j]
j - 1
ts.plot(na.omit(Hack[,j,with=FALSE]))
j <- j+1

series_names <- c("ptp_p",                       "ptp_f",                       "trend_weight",               
                  "idx",                         "unique_active_job_count",     "build_avg_time_to_fill",      "build_avg_salary",           
                  "build_n_postings",            "restaurant_avg_time_to_fill", "restaurant_avg_salary",       "restaurant_n_postings",      
                  "car_asi",                     "30_car_sai",                  "24_car_sai",                  "12_car_sai",                 
                  "air_asi")

needs_sa = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

take_logs =  c(0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
take_diffs = c(1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0)

hack_lib <- data.table(series_name = series_names, needs_SA = needs_sa, 
                       take_logs = take_logs, take_diffs = take_diffs)

write.csv(hack_lib, "C:\\Users\\seton\\Dropbox\\Quantagon\\inflation_nowcasting\\hack_lib.csv", row.names = FALSE)




