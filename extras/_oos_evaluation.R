# Out of sample evaluation -----------------------------------------------------
rm(list=ls()) # clear workspace

library(tidyverse)
library(tsbox)
library(data.table)

library(ggthemes)
long_names <- tribble(
  ~long , ~id,
  "CPI"  , "cpi",
  "CPI Core"  , "cpi_core",
  "PCE"  , "pce",
  "PCE Core"  , "pce_core")

long_models <- tribble(
  ~lmodel , ~model,
  "Cleveland Fed Nowcast"  , "cleveland fed NC",
  "Random Forest"  , "rf",
  "Random Forest, random node"  , "rfrn",
  "AR(1)-Benchmark", "a1")

setwd(here::here())
setwd("../")
main_dir<- getwd()

# --- Fcst eval function -------------------------------------------------------
# Diebold Mariano Test
dm.test <- function (e1, e2, h = 1, power = 2) {
  d <- c(abs(e1))^power - c(abs(e2))^power
  d.cov <- acf(d, na.action = na.omit, lag.max = h - 1,
               type = "covariance", plot = FALSE)$acf[, , 1]
  d.var <- sum(c(d.cov[1], 2 * d.cov[-1]))/length(d)
  dv <- d.var#max(1e-8,d.var)
  if(dv > 0)
    STATISTIC <- mean(d, na.rm = TRUE) / sqrt(dv)
  else if(h==1)
    stop("Variance of DM statistic is zero")
  else
  {
    warning("Variance is negative, using horizon h=1")
    return(dm.test(e1,e2,alternative,h=1,power))
  }
  n <- length(d)
  k <- ((n + 1 - 2*h + (h/n) * (h-1))/n)^(1/2)
  STATISTIC <- STATISTIC * k
  names(STATISTIC) <- "DM"
}

# Evaluation period ------------------------------------------------------------

start_date <- "2013-10-01"
# end_date <- "2022-01-01"
end_date <- "2019-12-31"
# 0) Load target data -----------------------------------------------------
DT <- fread(paste0(main_dir, "/_data/usa_data_sa.csv")) 
lib <- fread(paste0(main_dir, "/_data/library.csv"))
lib <- lib[country == "united states"]
DT[ , ref_date := as.Date(ref_date)]
DT[ , pub_date := as.Date(pub_date)]

meta <- tribble(
  ~series_name , ~id,
  "consumer price index cpi", "cpi",
  "core consumer prices", "cpi_core",
  "core pce price index", "pce_core",
  "personal consumption expenditures price index", "pce"
)

dta <- DT %>% 
  filter(series_name %in% c("consumer price index cpi",
                            "core consumer prices",
                            "core pce price index",
                            "personal consumption expenditures price index")
  ) %>% 
  left_join(meta, by = "series_name") %>% 
  select(pub_lag, id, time = ref_date, value) %>% 
  ts_span(start = "1990") 

final_mom <- dta %>% 
  select(-pub_lag) %>% 
  mutate(year = year(time), 
         month = month(time),
         time = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(id, ref_date = time, final = value) %>% 
  ts_pc()

final_yoy <- dta %>% 
  select(-pub_lag) %>% 
  mutate(year = year(time), 
         month = month(time),
         time = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(id, ref_date = time, final = value) %>% 
  ts_pcy()

# 1) Load univariate benchmarks ------------------------------------------------

load(paste0(main_dir, "/_results/benchmarks.Rdata"))
benchmarks <- out

# 2) Load clevaland fed nowcast ------------------------------------------------
load(paste0(main_dir, "/_results/cleveland_nc.Rdata"))
cleveland <- out

conv_names <- tribble(
  ~id , ~conv,
  "CPI.Inflation"  , "cpi",
  "Core.CPI.Inflation"  , "cpi_core",
  "PCE.Inflation"  , "pce",
  "Core.PCE.Inflation"  , "pce_core")

cleve_mom <- cleveland$mom %>% 
  left_join(conv_names, by = "id") %>% 
  select(id = conv, pub_date, ref_date, value)

cleve_yoy <- cleveland$yoy %>% 
  left_join(conv_names, by = "id") %>% 
  select(id = conv, pub_date, ref_date, value)

# 3) Load Random Forest nowcasts -----------------------------------------------
##  a) CPI ---------------------------------------------------------------------
load(paste0(main_dir, "/_results/cpi.rf_fit_part1.Rdata"))
cpi_rf.p1 <- cpi.rf_fit

load(paste0(main_dir, "/_results/cpi.rfrn_fit_part1.Rdata"))
cpi_rfrn.p1 <- cpi.rfrn_fit

load(paste0(main_dir, "/_results/cpi.rf_fit_part2.Rdata"))
cpi_rf.p2 <- cpi.rf_fit

load(paste0(main_dir, "/_results/cpi.rfrn_fit_part2.Rdata"))
cpi_rfrn.p2 <- cpi.rfrn_fit

# merge the files

cpi_oos_rf <- cpi_rf.p1 %>% 
  bind_rows(cpi_rf.p2 %>% filter(pub_date > "2019-10-19")) %>% 
  mutate(id = "cpi", 
         model = "rf")

cpi_oos_rfrn <- cpi_rfrn.p1 %>% 
  bind_rows(cpi_rfrn.p2 %>% filter(pub_date > "2019-10-19")) %>% 
  mutate(id = "cpi", 
         model = "rfrn")

# Combine results 
cpi.rf_oos <- bind_rows(cpi_oos_rf, cpi_oos_rfrn) %>% 
  mutate(year = year(ref_date), 
         month = month(ref_date),
         ref_date = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(model, id, pub_date, ref_date, value) 

##  b) CPI CORE ----------------------------------------------------------------
load(paste0(main_dir, "/_results/cpi_core.rf_fit_part1.Rdata"))
cpi_core.rf_fit.p1 <- cpi_core.rf_fit

load(paste0(main_dir, "/_results/cpi_core.rfrn_fit_part1.Rdata"))
cpi_core.rfrn_fit.p1 <- cpi_core.rfrn_fit

load(paste0(main_dir, "/_results/cpi_core.rf_fit_part2.Rdata"))
cpi_core.rf_fit.p2 <- cpi_core.rf_fit

load(paste0(main_dir, "/_results/cpi_core.rfrn_fit_part2.Rdata"))
cpi_core.rfrn_fit.p2 <- cpi_core.rfrn_fit

load(paste0(main_dir, "/_results/cpi_core.rf_fit_part3.Rdata"))
cpi_core.rf_fit.p3 <- cpi_core.rf_fit

load(paste0(main_dir, "/_results/cpi_core.rfrn_fit_part3.Rdata"))
cpi_core.rfrn_fit.p3 <- cpi_core.rfrn_fit

cpi_core.rf_oos <- bind_rows(cpi_core.rf_fit.p1, 
                             cpi_core.rf_fit.p2,
                             cpi_core.rf_fit.p3) %>% 
  mutate(year = year(ref_date), 
         month = month(ref_date),
         ref_date = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(pub_date, ref_date, value) %>% 
  mutate(id = "cpi_core", 
         model = "rf")

cpi_core.rfrn_oos <- bind_rows(cpi_core.rfrn_fit.p1, 
                             cpi_core.rfrn_fit.p2,
                             cpi_core.rfrn_fit.p3) %>% 
  mutate(year = year(ref_date), 
         month = month(ref_date),
         ref_date = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(pub_date, ref_date, value) %>% 
  mutate(id = "cpi_core", 
         model = "rfrn")

# Combine RF NC

cpi_core.rf_oos <- bind_rows(cpi_core.rf_oos, 
            cpi_core.rfrn_oos)

##  c) PCE ----------------------------------------------------------------
load(paste0(main_dir, "/_results/pce.rf_fit_part1.Rdata"))
pce.rf_fit.p1 <- pce.rf_fit

load(paste0(main_dir, "/_results/pce.rfrn_fit_part1.Rdata"))
pce.rfrn_fit.p1 <- pce.rfrn_fit

load(paste0(main_dir, "/_results/pce.rf_fit_part2.Rdata"))
pce.rf_fit.p2 <- pce.rf_fit

load(paste0(main_dir, "/_results/pce.rfrn_fit_part2.Rdata"))
pce.rfrn_fit.p2 <- pce.rfrn_fit

load(paste0(main_dir, "/_results/pce.rf_fit_part3.Rdata"))
pce.rf_fit.p3 <- pce.rf_fit

load(paste0(main_dir, "/_results/pce.rfrn_fit_part3.Rdata"))
pce.rfrn_fit.p3 <- pce.rfrn_fit

pce.rf_oos <- bind_rows(pce.rf_fit.p1, 
                             pce.rf_fit.p2,
                             pce.rf_fit.p3) %>% 
  mutate(year = year(ref_date), 
         month = month(ref_date),
         ref_date = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(pub_date, ref_date, value) %>% 
  mutate(id = "pce", 
         model = "rf")

pce.rfrn_oos <- bind_rows(pce.rfrn_fit.p1, 
                               pce.rfrn_fit.p2,
                               pce.rfrn_fit.p3) %>% 
  mutate(year = year(ref_date), 
         month = month(ref_date),
         ref_date = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(pub_date, ref_date, value) %>% 
  mutate(id = "pce", 
         model = "rfrn")

# Combine RF NC

pce.rf_oos <- bind_rows(pce.rf_oos, 
                             pce.rfrn_oos)

##  d) PCE CORE ----------------------------------------------------------------
load(paste0(main_dir, "/_results/pce_core.rf_fit_part1.Rdata"))
pce_core.rf_fit.p1 <- pce_core.rf_fit

load(paste0(main_dir, "/_results/pce_core.rfrn_fit_part1.Rdata"))
pce_core.rfrn_fit.p1 <- pce_core.rfrn_fit

load(paste0(main_dir, "/_results/pce_core.rf_fit_part2.Rdata"))
pce_core.rf_fit.p2 <- pce_core.rf_fit

load(paste0(main_dir, "/_results/pce_core.rfrn_fit_part2.Rdata"))
pce_core.rfrn_fit.p2 <- pce_core.rfrn_fit

load(paste0(main_dir, "/_results/pce_core.rf_fit_part3.Rdata"))
pce_core.rf_fit.p3 <- pce_core.rf_fit

load(paste0(main_dir, "/_results/pce_core.rfrn_fit_part3.Rdata"))
pce_core.rfrn_fit.p3 <- pce_core.rfrn_fit

pce_core.rf_oos <- bind_rows(pce_core.rf_fit.p1, 
                        pce_core.rf_fit.p2,
                        pce_core.rf_fit.p3) %>% 
  mutate(year = year(ref_date), 
         month = month(ref_date),
         ref_date = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(pub_date, ref_date, value) %>% 
  mutate(id = "pce_core", 
         model = "rf")

pce_core.rfrn_oos <- bind_rows(pce_core.rfrn_fit.p1, 
                          pce_core.rfrn_fit.p2,
                          pce_core.rfrn_fit.p3) %>% 
  mutate(year = year(ref_date), 
         month = month(ref_date),
         ref_date = as.Date(paste0(year, "-", month, "-01"))) %>% 
  select(pub_date, ref_date, value) %>% 
  mutate(id = "pce_core", 
         model = "rfrn")

# Combine RF NC

pce_core.rf_oos <- bind_rows(pce_core.rf_oos, 
                        pce_core.rfrn_oos)


##  e) Bring RF NC together -----------------------------------------------------

rf_oos <-bind_rows(cpi.rf_oos,
                   cpi_core.rf_oos,
                   pce.rf_oos,
                   pce_core.rf_oos)

# 4) Compute errors ------------------------------------------------------------
##  a) Benchmark ---------------------------------------------------------------
# MoM
benchmarks_errors.mom <- benchmarks$mom %>% 
  separate(expr, into = c("id", "model"), sep = "\\.") %>% 
  left_join(final_mom, by = c("id","ref_date")) %>% 
  filter(ref_date > start_date) %>% 
  filter(ref_date < end_date) %>% 
  group_by(id, model, ref_date) %>% 
  arrange(desc(pub_date)) %>% 
  mutate(h = seq(n())) %>% 
  ungroup() %>%
  mutate(error = final - value) %>%
  filter(!is.na(error)) %>%
  mutate(h = as.numeric(h)) %>% 
  # mutate(h = if_else(id == "pce" | id == "pce_core", h-4,h-2)) %>% 
  filter(h < 48)
  
# YoY
benchmarks_errors.yoy <- benchmarks$yoy %>% 
  separate(expr, into = c("id", "model"), sep = "\\.") %>% 
  left_join(final_yoy, by = c("id","ref_date")) %>% 
  filter(ref_date > start_date) %>% 
  filter(ref_date < end_date) %>% 
  group_by(id, model, ref_date) %>% 
  arrange(desc(pub_date)) %>% 
  mutate(h = seq(n())) %>% 
  ungroup() %>%
  mutate(error = final - value) %>%
  filter(!is.na(error)) %>%
  mutate(h = as.numeric(h)) %>% 
  # mutate(h = if_else(id == "pce" | id == "pce_core", h-4,h-2)) %>% 
  filter(h < 48)

##  b) Cleveland Fed ------------------------------------------------------------
# Cleveland NC is daily, we select the same days as in benchmarks
# MoM 
bench_pub_dates <- benchmarks_errors.mom %>% distinct(pub_date) 
bench_pub_dates <- bench_pub_dates$pub_date - lubridate::days(1) # want fridays not saturday

cleve_errors.mom <- cleve_mom %>% 
  left_join(final_mom, by = c("id","ref_date")) %>% 
  filter(ref_date > start_date) %>% 
  filter(ref_date < end_date) %>% 
  filter(pub_date %in% bench_pub_dates) %>% 
  group_by(id, ref_date) %>% 
  arrange(desc(pub_date)) %>% 
  mutate(h = seq(n())) %>% 
  ungroup() %>%
  mutate(error = final - value) %>%
  filter(!is.na(error)) %>%
  mutate(h = as.numeric(h))  %>% 
  # mutate(h = case_when(
  #   id == "pce" ~ h-4,
  #   id == "pce_core" ~ h-4,
  #   id == "cpi" ~ h-2,
  #   id == "cpi_core" ~ h-2
  # )) %>% 
  # mutate(h = if_else(ref_date >= "2018-12-01" & ref_date <= "2019-02-01", 
  #                    case_when(
  #                      id == "pce" ~ h-3,
  #                      id == "pce_core" ~ h-3
  #                    ), 
  #                    h
  # )
  # ) %>% 
  mutate(pub_date = pub_date + lubridate::days(1))

cleve_errors.yoy <- cleve_yoy %>% 
  left_join(final_yoy, by = c("id","ref_date")) %>% 
  filter(ref_date > start_date) %>% 
  filter(ref_date < end_date) %>% 
  filter(pub_date %in% bench_pub_dates) %>% 
  group_by(id, ref_date) %>% 
  arrange(desc(pub_date)) %>% 
  mutate(h = seq(n())) %>% 
  ungroup() %>%
  mutate(error = final - value) %>%
  filter(!is.na(error)) %>%
  mutate(h = as.numeric(h))  %>% 
  # mutate(h = case_when(
  #   id == "pce" ~ h-4,
  #   id == "pce_core" ~ h-4,
  #   id == "cpi" ~ h-2,
  #   id == "cpi_core" ~ h-2
  # )) %>% 
  # mutate(h = if_else(ref_date >= "2018-12-01" & ref_date <= "2019-02-01", 
  #                    case_when(
  #                      id == "pce" ~ h-3,
  #                      id == "pce_core" ~ h-3
  #                      ), 
  #                    h
  #                    )
  #        ) %>% 
  mutate(pub_date = pub_date - lubridate::days(2))



##  c) RF Models ---------------------------------------------------------------

rf_errors <- rf_oos %>% 
  left_join(final_mom, by = c("id","ref_date")) %>% 
  filter(ref_date > start_date) %>% 
  filter(ref_date < end_date) %>% 
  group_by(model, id, ref_date) %>% 
  arrange(desc(pub_date)) %>% 
  mutate(h = seq(n())) %>% 
  ungroup() %>%
  mutate(error = final - value) %>%
  filter(!is.na(error)) %>%
  mutate(h = as.numeric(h)) # %>% 
  # mutate(h = case_when(
  #   id == "pce" ~ h-4,
  #   id == "pce_core" ~ h-4,
  #   id == "cpi" ~ h-2,
  #   id == "cpi_core" ~ h-2
  # ))

# 5) Visual inspection ------------------------------------------------------------------------
# Benchmarks
benchmarks_errors.mom %>%
  group_by(id, model, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(id, model, h, rmse) %>%
  ggplot() +
  geom_line(aes(y = rmse,
                x = h,
                color = model,
                linetype = id),
  )
# Plot data
benchmarks_errors.yoy %>%
  group_by(id, model, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(id, model, h, rmse) %>%
  ggplot() +
  geom_line(aes(y = rmse,
                x = h,
                color = model,
                linetype = id),
  )


# Plot data
cleve_errors.yoy %>%
  group_by(id, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(id, h, rmse) %>%
  filter(h >= -3 & h <= 4) %>% 
  ggplot() +
  geom_line(aes(y = rmse,
                x = h,
                linetype = id),
  )


# Plot data
cleve_errors.mom %>%
  group_by(id, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(id, h, rmse) %>%
  filter(h >= -3 & h <= 4) %>% 
  ggplot() +
  geom_line(aes(y = rmse,
                x = h,
                linetype = id),
  )

# Plot data
rf_errors %>%
  group_by(model, id, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(model, id, h, rmse) %>%
  filter(h >= -3 & h <= 3) %>%
  ggplot() +
  geom_line(aes(y = rmse,
                x = h,
                linetype = id, 
                color = model),
  )

# 6) Obtains RMSE --------------------------------------------------------------

# AR1-benchmark
bench_rmse <- benchmarks_errors.mom %>% 
  group_by(id, model, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(id, model, h, rmse) %>%
  filter(model == "a1") 
  
# Cleveland MoM
cleve_rmse <- cleve_errors.mom %>% 
  group_by(id, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(id, h, rmse) %>% 
  filter(h <= 5) %>% 
  mutate(model = "cleveland fed NC")

# RF Models
rf_rmse  <- rf_errors %>% 
  group_by(id, model, h) %>%
  summarize(mse = sum(error^2),
            rmse = sqrt(sum(error^2)),
            mae = mean(abs(error))) %>%
  ungroup() %>%
  select(id, model, h, rmse) %>% 
  filter(h <= 4) 

# 7) Main Figure ---------------------------------------------------------------
# GRAPH 1: BRIDGE MODEL COMPARISON

fig1 <- bind_rows(cleve_rmse, rf_rmse) %>% 
  left_join(bench_rmse, by = c("id", "h")) %>% 
  mutate(rrmse = rmse.x/rmse.y) %>% 
  # filter(h >= -3 & h <= 5) %>%
  select(id, h, rrmse, model = model.x) %>% 
  left_join(long_names, by = "id") %>% 
  left_join(long_models, by = "model") %>% 
  ggplot() +
  geom_line(aes(y = rrmse,
                x = h,
                color = lmodel),
  ) + 
  facet_wrap(vars(long), scales = "fixed") +
  geom_hline(yintercept=1 , colour="#999999", size = .5, alpha = 0.5) +
  geom_vline(xintercept=c(0) , colour="#999999", size = .5, alpha = 0.5) +
  geom_text(aes(x=0, label="Backcast", y=1.45), hjust = 1.5, colour="#000000", family = "serif")+
  geom_text(aes(x=3, label="Nowcast", y=1.45), hjust = 1.5, colour="#000000", family = "serif")+
  theme_tufte(base_size = 12, base_family = "serif", ticks = TRUE) +
  scale_color_manual(values = c("#999999","#000000", "#7F7F7F")) +
  labs(x = "Horizon in weeks",
       y = "Relative RMSE")+
  theme(axis.text = element_text(size=10))+
  scale_y_continuous(
    expand = expand_scale(0.005),
    limits = c(0, 1.5),
    breaks = seq(0, 1.5, by = 0.25)) +
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  theme(legend.text=element_text(size=10),
        legend.position = "bottom",
        legend.title = element_blank())
fig1
ggsave(paste0(main_dir,"/_results/","rrmse.png"), 
       width = 8, height = 4.5, dpi = 300, units = "in", device = "png")

dev.off()

# 8) Main Table and Diebold Mariano Test ---------------------------------------
## a) Table with relative RMSE -------------------------------------------------
tab_bench <- bench_rmse %>% 
  filter(model == "a1") %>% 
  # filter(case_when(
  #   id == "pce" ~ h >= -3 & h <= 0,
  #   id == "pce_core" ~ h >= -3 & h <= 0,
  #   id == "cpi" ~ h >= -3 & h <= 2,
  #   id == "cpi_core" ~ h >= -3 & h <= 2)) %>% 
  filter(h <= 4) %>% 
  left_join(long_names, by = "id") %>% 
  left_join(long_models, by = "model") %>% 
  select(-id, -model) %>% 
  group_by(long) %>% 
  arrange(h) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "h", values_from = "rmse") 
  
tab <- bind_rows(cleve_rmse, rf_rmse) %>% 
  left_join(bench_rmse, by = c("id", "h")) %>% 
  mutate(rrmse = rmse.x/rmse.y) %>% 
  # filter(case_when(
  #   id == "pce" ~ h >= -3 & h <= 0,
  #   id == "pce_core" ~ h >= -3 & h <= 0,
  #   id == "cpi" ~ h >= -3 & h <= 2,
  #   id == "cpi_core" ~ h >= -3 & h <= 2)) %>%
  filter(h <= 4) %>% 
  select(id, h, rrmse, model = model.x) %>% 
  left_join(long_names, by = "id") %>% 
  left_join(long_models, by = "model") %>%
  select(-id, -model) %>% 
  group_by(long) %>% 
  arrange(h) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "h", values_from = "rrmse") %>% 
  bind_rows(tab_bench) %>% 
  mutate_if(is.numeric, round, 3) %>% 
  arrange(long, lmodel) %>% 
  rename(`model` = "lmodel",
         `target` = "long")
  
  

# paste the code in the latex document and adapt it
library(kableExtra)
tab %>%
  kable(format = "latex",
        booktabs = TRUE,
        linesep = "") %>%
  kable_styling(font_size = 9, full_width = F) 


## b) Diebold Mariano Test -----------------------------------------------------

# H0: Nowcast is better than AR(1) benchmark

b.e <- benchmarks_errors.mom %>%
  filter(h > 0) %>%
  filter(h < 5) %>% 
  filter(model == "a1")

cl.e <- cleve_errors.mom %>%
  filter(h > 0) %>%
  filter(h < 5) %>% 
  mutate(model = "cleveland fed NC")

rf.e <- rf_errors %>%
  filter(h > 0) %>%
  filter(h < 5)

mod.e <- bind_rows(cl.e, rf.e)

# Q 1: Are nowcasts better than AR(1)-Benchmark? 

# Create a function that iterate over all time horizons

dm.test_h <-function(x, y){
  
  res.1 <- list()
  
  for (i in 1:4) {
    e.1 <- b.e %>%
      filter(h == i & id == y, model == "a1") %>%
      select(error) %>%
      deframe()
    
    e.2 <- mod.e %>%
      filter(model == x & id == y & h == i) %>%
      select(error) %>%
      deframe()
    
    res.1[[i]] <- tibble("mod.1" = "a1",
                         "mod.2" = x,
                         "target" = y,
                         "horizon" = i,
                         test_value = round(forecast::dm.test(e.1,
                                                              e.2,
                                                              alternative = "greater" ,
                                                              h = i, power = 2)$statistic, 2),
                         p_value  = round(forecast::dm.test(e.1,
                                                            e.2,
                                                            alternative = "greater",
                                                            h = i, power = 2)$p.value, 3))}
  test_res1 <- res.1 %>% bind_rows() 
  
}

# use lapply to loop over all models

models <- mod.e %>% 
  distinct(model) %>% 
  deframe() %>% 
  as.list()

targets <- mod.e %>% 
  distinct(id) %>% 
  deframe() %>% 
  as.list()

results.cpi <- lapply(models, dm.test_h, y = "cpi") %>% bind_rows()
results.cpi_core <- lapply(models, dm.test_h, y = "cpi_core") %>% bind_rows()
results.pce <- lapply(models, dm.test_h, y = "pce") %>% bind_rows()
results.pce_core <- lapply(models, dm.test_h, y = "pce_core") %>% bind_rows()

p.value_fct <- function (p.value) {
  unclass(symnum(p.value, corr = FALSE, na = FALSE, 
                 cutpoints = c(0, 0.01, 0.05, 0.1, 1), 
                 symbols = c("***", "**", "*", " ")))}


tab_res1 <- results.cpi %>% 
  bind_rows(results.cpi_core, results.pce, results.pce_core) %>% 
  mutate(stars = p.value_fct(p_value)) %>% 
  unite(value, c(test_value, stars), sep = "") %>% 
  select(-p_value) %>% 
  pivot_wider(names_from = "horizon", values_from = "value") %>% 
  rename(model = mod.2, 
         id = target) %>% 
  select(-mod.1) %>% 
  left_join(long_names, by = "id") %>% 
  left_join(long_models, by = "model") %>% 
  arrange(long, lmodel) %>% 
  select(-model, -id) %>% 
  rename(`model` = "lmodel",
         `target` = "long") %>% 
  select(target, model, everything())

tab_res1 %>% 
  kable(format = "latex",
        booktabs = TRUE, 
        linesep = "") %>%
  kable_styling(font_size = 9, full_width = F) %>% 
  add_header_above(c(" " = 1, "Horizon" = 5))


# Q 2: Is RF better than Cleveland Fed? 

# Create a function that iterate over all time horizons

dm.test_h <-function(x, y){
  
  res.1 <- list()
  
  for (i in 1:4) {
    e.1 <- mod.e %>%
      filter(h == i & id == y, model == "cleveland fed NC") %>%
      select(error) %>%
      deframe()
    
    e.2 <- mod.e %>%
      filter(model == x & id == y & h == i) %>%
      select(error) %>%
      deframe()
    
    res.1[[i]] <- tibble("mod.1" = "cleveland fed NC",
                         "mod.2" = x,
                         "target" = y,
                         "horizon" = i,
                         test_value = round(forecast::dm.test(e.1,
                                                              e.2,
                                                              alternative = "greater" ,
                                                              h = i, power = 2)$statistic, 2),
                         p_value  = round(forecast::dm.test(e.1,
                                                            e.2,
                                                            alternative = "greater",
                                                            h = i, power = 2)$p.value, 3))}
  test_res1 <- res.1 %>% bind_rows() 
  
}

# use lapply to loop over all models

models <- mod.e %>% 
  distinct(model) %>% 
  deframe() %>% 
  as.list()

results.cpi <- lapply(models[2:3], dm.test_h, y = "cpi") %>% bind_rows()
results.cpi_core <- lapply(models[2:3], dm.test_h, y = "cpi_core") %>% bind_rows()
results.pce <- lapply(models[2:3], dm.test_h, y = "pce") %>% bind_rows()
results.pce_core <- lapply(models[2:3], dm.test_h, y = "pce_core") %>% bind_rows()

tab_res2 <- results.cpi %>% 
  bind_rows(results.cpi_core, results.pce, results.pce_core) %>% 
  mutate(stars = p.value_fct(p_value)) %>% 
  unite(value, c(test_value, stars), sep = "") %>% 
  select(-p_value) %>% 
  pivot_wider(names_from = "horizon", values_from = "value") %>% 
  rename(model = mod.2, 
         id = target) %>% 
  select(-mod.1) %>% 
  left_join(long_names, by = "id") %>% 
  left_join(long_models, by = "model") %>% 
  arrange(long, lmodel) %>% 
  select(-model, -id) %>% 
  rename(`model` = "lmodel",
         `target` = "long") %>% 
  select(target, model, everything())

tab_res2 %>% 
  kable(format = "latex",
        booktabs = TRUE, 
        linesep = "") %>%
  kable_styling(font_size = 9, full_width = F) %>% 
  add_header_above(c(" " = 1, "Horizon" = 5))
