library(mlbench)
library(caret)
library(BART)
library(RegTree)

data(BostonHousing)
dim(BostonHousing)
y <- BostonHousing$medv
df <- within(BostonHousing, rm(medv))
df <- sapply(df,as.numeric)

set.seed(42) 
test_inds = createDataPartition(y = 1:length(y), p = 0.2, list = F) 
df_test = df[test_inds, ] 
y_test = y[test_inds] 
df_train = df[-test_inds, ]
y_train = y[-test_inds]
# fit model

model <- wbart(x.train = df_train, y.train = y_train, x.test = df_test)
ts.plot(cbind(y_test,model$yhat.test.mean), col=c("black", "red"))

Btree <- reg_forest(y_train, df_train, max_nodes = 41, type = "bayes")

Bfit <- Btree$in_sample_fit
ts.plot(cbind(y_train,Bfit), col=c("black", "red"))
