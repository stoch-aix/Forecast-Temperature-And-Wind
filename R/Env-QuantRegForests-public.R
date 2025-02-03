
#
# Prediction using quantile regression forests
#
# code based on a template from Gemini 2.0 flash
# 
# (C) A. Steland, Institute of Statistics and AI Center, RWTH Aachen University
# published under the Academic Research Software Licence (ARSL) 1.0
#
# 
 

library(forecast)
library(Hmisc)
library(dplyr)
library(fastDummies)
library(scoringutils)

# Load libraries
library(quantregForest)
library(caret)
library(dplyr)

source("helpers.R")

#### Get data and divide into training and test

load("ICData.rda")
dat=station_int[station_int[,4]==12,]
dat=data.frame(dat)

dat = add_nonlinearities(dat)

dat2 <- fastDummies::dummy_cols(dat,select_columns = "month",remove_first_dummy = TRUE)
dat2 <-dat2 %>% mutate(id=row_number())
train=dat2 %>% filter((month==9 & day<22 & year==2021)| (month<9 & year==2021) | (year<2021))
test= anti_join(dat2, train,by='id')

### Create the regression variables

# only using the past up to lag 3: MSE = 0.4754
regressors.wind = c( "windspeedmph", "wind.l1", "wind.l2", "wind.l3" )

regressors.wind = c( "humidity", "windspeedmph", "solarradiation", "tempf", "wind.l1", "wind.l2", "wind.l3" )

# gives MSE: 49.86
regressors.temp = c("humidity", "solarradiation", "temp.log", "tempf",  "tempf2", "temp.l1", "temp.l2" )

regressors.temp = c("humidity", "solarradiation", "tempf",  "windspeedmph", "temp.l1", "temp.l2", "temp.l3" )


# create learning samples

xreg.wind.train=data.matrix(train[,c(regressors.wind)])
xreg.temp.train=data.matrix(train[,regressors.temp])

# create test sample

xreg.wind.test=test[,regressors.wind]
xreg.temp.test=test[,regressors.temp]


# temperature  (omit first obs with NAs due to lagged values)
X.train = xreg.temp.train[4:nrow(xreg.temp.train),]
X.test = xreg.temp.test[1:100,]
y.train = train$tmpf.res[4:nrow(xreg.temp.train)]
y.test = test$tmpf.res[1:100] # there is a NA at the end... 


# Quantile Regression Forests

# for (approx.) reproducible results, one needs to restart R 
set.seed(123456)

# Train the Quantile Regression Forest model
qrf_temp_model <- quantregForest(
  x = X.train,
  y = y.train, nthreads = 10,
  ntree = 1000  # Number of trees; increase for better performance
)

# Predict the mean temperature, predict whole test sample using training fit
predictions <- predict(qrf_temp_model, newdata = X.test, what = c(0.025, 0.5, 0.975))
# put to format pred, q_lower, q_upper
pred.temp = cbind( predictions[,2], predictions[,1], predictions[,3] )

# sequential predictions

pred.temp=matrix(nrow=(nrow(test)-1),ncol=3)

qrf.temp.model2 = quantregForest(
  x = X.train,
  y = y.train, nthreads = 10,
  ntree = 1000  # Number of trees; increase for better performance
)

predictions <- predict(qrf.temp.model2, newdata = X.test[1,], what = c(0.025, 0.5, 0.975))
pred.temp[1,] = c( predictions[,2], predictions[,1], predictions[,3] )


### Forecast
for (i in 2:nrow(pred.wind)){
  qrf.temp.model2 = quantregForest( y=c(y.train,y.test[1:(i-1)]), 
    x=as.matrix(rbind(X.train, X.test[1:(i-1),])),ntree=1000, nthreads=10 )
  predictions <- predict(qrf.temp.model2, newdata = X.test[i,], what = c(0.025, 0.5, 0.975))
  pred.temp[i,] = c( predictions[,2], predictions[,1], predictions[,3] )
}
dput( pred.temp, "qrf-wind-prediction-intervals.dat" )


# wind speed

# wind 
X.train = xreg.wind.train[4:nrow(xreg.wind.train),]
X.test = xreg.wind.test[1:100,]
y.train = train$wind.res[4:length(train$wind.res)]
y.test = test$wind.res[1:100] # there is a NA at the end... 

# Train the Quantile Regression Forest model
# short: 
# fit = quantregForest( x = X.train, y = y.train, ntree=1000, what = c(0.5) )
# ypred = predict( fit, newdata = , what = )
qrf_wind_model <- quantregForest(
  x = X.train,
  y = y.train,
  ntree = 1000  # Number of trees; increase for better performance
)

# Predict the mean temperature
predictions <- predict(qrf_wind_model, newdata = X.test, what = c(0.025, 0.5, 0.975))
# put to format pred, q_lower, q_upper
pred.wind = cbind( predictions[,2], predictions[,1], predictions[,3] )

# sequential predictions
pred.wind=matrix(nrow=(nrow(test)-1),ncol=3)

qrf.wind.model2 = quantregForest(
  x = X.train,
  y = y.train, nthreads = 10,
  ntree = 1000  # Number of trees; increase for better performance
)

predictions <- predict(qrf.wind.model2, newdata = X.test[1,], what = c(0.025, 0.5, 0.975))
pred.wind[1,] = c( predictions[,2], predictions[,1], predictions[,3] )


### Forecast
for (i in 2:nrow(pred.wind)){
  qrf.wind.model2 = quantregForest( y=c(y.train,y.test[1:(i-1)]), 
                                    x=as.matrix(rbind(X.train, X.test[1:(i-1),])),ntree=1000, nthreads=10 )
  predictions <- predict(qrf.wind.model2, newdata = X.test[i,], what = c(0.025, 0.5, 0.975))
  pred.wind[i,] = c( predictions[,2], predictions[,1], predictions[,3] )
}
dput( pred.wind, "qrf-wind-prediction-intervals.dat" )

# Analyze results

### Actual values to be predicted (one-day ahead values)
actual.temp=test$tmpf.res[1:(nrow(test)-1)]
actual.wind=test$wind.res[1:(nrow(test)-1)]

# the results

tab = generate_table( actual.temp, actual.wind, pred.temp, pred.wind ) 
tab

CI_MSE( actual.temp-pred.temp[,1] )

dput( actual.temp-pred.temp[,1], "qrf-temp-residuals.dat" )
dput( actual.wind-pred.wind[,1], "qrf-wind-residuals.dat" )

cat("Test temperature models...\n")
Test_Models_MSE( dget("qrf-temp-residuals.dat"), dget("arima-temp-residuals.dat"), alpha=0.05 )
cat("Test wind models...\n")
Test_Models_MSE( dget("qrf-wind-residuals.dat"), dget("arima-wind-residuals.dat"), alpha=0.05 )



