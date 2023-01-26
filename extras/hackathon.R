library(RegTree)
library(dateutils)
library(devtools)
library(data.table)
# load_all()

dt <- fread("/tmp/tickerdata.csv")
dt[ , as_of_date := as.Date(as_of_date)]

train <- dt[dt$as_of_date < as.Date("2021-01-01")]
test <- dt[dt$as_of_date >= as.Date("2021-01-01")]

y <- train$IDXX
X <- as.matrix(train[,-c(1,2),with=FALSE])
ytest <- test$IDXX
Xtest <- as.matrix(test[,-c(1,2),with=FALSE])


Btree <- reg_forest(y, X, max_nodes = 41, draws = 5000, type = "bayes", weight_pow = 2)
Rforest <- reg_forest(y, X, max_nodes = 81)

yf <- Btree$in_sample_fit
dt_insamp <- cbind(y, 3*(yf-mean(yf)) + mean(yf))
pretty_plot(dt_insamp)

yf <- Rforest$out_of_sample
rf_insamp <- cbind(y, yf)
pretty_plot(rf_insamp)

tst <- StdFitField(Xtest, Btree$Trees)
yt <- tst[[1]]
dt_test <- cbind(ytest, yt)
pretty_plot(dt_test)