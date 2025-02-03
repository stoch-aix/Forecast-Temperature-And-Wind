#
# Re-analysis 
#


#### Get data and divide into training and test

load("ICData.rda")

dat=station_int[station_int[,4]==12,]
library(forecast)
library(Hmisc)
library(dplyr)
library(lightgbm)
library(forecast)
library(fastDummies)
library(sandwich)
library(scoringutils)

source("helpers.R")

dat=data.frame(dat)

dat$wind.res=Lag(dat$windspeedmph,-1)
dat$tmpf.res=Lag(dat$tempf,-1)

dat2 <- fastDummies::dummy_cols(dat,select_columns = "month",remove_first_dummy = TRUE)
dat2 <-dat2 %>% mutate(id=row_number())
train=dat2 %>% filter((month==9 & day<22 & year==2021)| (month<9 & year==2021) | (year<2021))
test = anti_join(dat2, train,by='id')

### set up regression models with ARIMA errors as in Bonas et al. (2024)

### Create the regression variables
xreg.wind.train=data.matrix(train[,c(5:7,12:22)])
xreg.temp.train=data.matrix(train[,c(6,12:22)])

xreg.wind.test=test[,c(5:7,12:22)]
xreg.temp.test=test[,c(6,12:22)]

###ARIMAX models

wind.model=auto.arima(train$wind.res,xreg = xreg.wind.train, seasonal = FALSE)
temp.model=auto.arima(train$tmpf.res,xreg = xreg.temp.train, seasonal = FALSE)

pred.wind=matrix(nrow=(nrow(test)-1),ncol=3)
pred.temp=matrix(nrow=(nrow(test)-1),ncol=3)

wind.model2=Arima(train$wind.res,xreg=as.matrix(xreg.wind.train),model=wind.model)
wind.1=data.frame(forecast(wind.model2,xreg=as.matrix(xreg.wind.test[1,]),h=1))
pred.wind[1,]=c(wind.1$Point.Forecast,wind.1$Lo.95,wind.1$Hi.95)

temp.model2=Arima(train$tmpf.res,xreg=as.matrix(xreg.temp.train),model=temp.model)
temp.1=data.frame(forecast(temp.model2,xreg=as.matrix(xreg.temp.test[1,]),h=1))
pred.temp[1,]=c(temp.1$Point.Forecast,temp.1$Lo.95,temp.1$Hi.95)

### Forecast
for (i in 2:nrow(pred.wind)){
  wind.model2=Arima(c(train$wind.res,test$wind.res[1:(i-1)]),xreg=as.matrix(rbind(xreg.wind.train,xreg.wind.test[1:(i-1),])),model=wind.model)
  wind.1=data.frame(forecast(wind.model2,xreg=as.matrix(xreg.wind.test[i,]),h=1))
  pred.wind[i,]=c(wind.1$Point.Forecast,wind.1$Lo.95,wind.1$Hi.95)
  
  temp.model2=Arima(c(train$tmpf.res,test$tmpf.res[1:(i-1)]),xreg=as.matrix(rbind(xreg.temp.train,xreg.temp.test[1:(i-1),])),model=temp.model)
  temp.1=data.frame(forecast(temp.model2,xreg=as.matrix(xreg.temp.test[i,]),h=1))
  pred.temp[i,]=c(temp.1$Point.Forecast,temp.1$Lo.95,temp.1$Hi.95)
  
}

### Actual values
actual.temp=test$tmpf.res[1:(nrow(test)-1)]
actual.wind=test$wind.res[1:(nrow(test)-1)]

### Replace predicted negative values with 0
pred.wind[pred.wind<0]<-0

### Actual values
actual.temp=test$tmpf.res[1:(nrow(test)-1)]
actual.wind=test$wind.res[1:(nrow(test)-1)]


# generate table with results

tab = generate_table( actual.temp, actual.wind, pred.temp, pred.wind ) 
tab


# test residuals for stationarity
library(tseries)
kpss.test( temp.model$residuals, null=c("Level") )
kpss.test( wind.model$residuals, null=c("Level") )

# ACF plot based on model residuals

par( mfrow = c(2,2) )
acf( temp.model$residuals, main="Temperature" )
acf( wind.model$residuals, main="Wind" ) 
pacf( temp.model$residuals, main=""  )
pacf( wind.model$residuals, main="" ) 

# ACF plot based on local demeaning 

wi = train$wind.res
wi.detrended = numeric( length(wi ) )
te = train$tempf
te.detrended = numeric( length(te) )
for ( i in (1:length(wi)) ) {
  h = c( 3, i, length(wi)-i+1 )
  h = min( h) 
  wi.detrended[i] = wi[i] - mean( wi[(i-h):i] )
  te.detrended[i] = te[i] - mean( te[(i-h):i] )
}
par( mfrow = c(2,2) )
acf( te.detrended, main="Temperature" )
acf( wi.detrended, main="Wind" ) 


# Extract parameters of the estimated model

ar1 = temp.model$coef[1]
b0 = temp.model$coef[2]
beta = temp.model$coef[3:length(fc$model$coef)]

# In-sample predictions 

yp = numeric(nrow(train))
eps = train$tempf - b0 - as.matrix(xreg.temp.train) %*% beta 
yy = c( train$tempf ) 
yp[1] = train$tempf[1]
for ( i in 2:nrow(train) ) {
  e = train$tempf[i] - b0 - crossprod( beta, as.numeric(xreg.temp.train[i,]) )
  #yp[i] = b0 + ar1 * eps[i] + crossprod( beta, as.numeric(xreg.temp.train[i,]) )
  yp[i] = b0 + crossprod( beta, as.numeric(xreg.temp.train[i,]) ) + ar1 * e
}
par( mfrow=c(1,1) )
plot( train$tempf, type="l", ylim=c(0,120) )
lines( yp, col="blue") 

# confidence interval for in-sample MSEP (training sample)
CI_MSE( train$tmpf.res - yp )


#
# Out-of-sample predictions based on trained model
# 
# These predictions slightly differ from the above sequential forecast in pred.temp,
# which are based on the sequentially refitted models. Contrary, here one fixes
# the model trained in the learning sample and applies it to the test sample.
# 

yp = numeric(nrow(test))
eps = test$tempf - b0 - as.matrix(xreg.temp.test) %*% beta 
yp[1] = test$tempf[1]
for ( i in 2:nrow(test) ) {
  e = test$tempf[i] - b0 - crossprod( beta, as.numeric(xreg.temp.test[i,]) )
  yp[i] = b0 + crossprod( beta, as.numeric(xreg.temp.test[i,]) ) + ar1 * e 
}
plot( test$tempf, type="l",  )
lines( yp, col="blue") 
# comparison: lines( pred.temp[,1], col="red")

# confidence interval for out-of-sample MSEP (test sample)
CI_MSE( test$tmpf.res[1:100] - yp[1:100] )

# save the sequential predictions...
#dput( actual.temp-pred.temp[,1], "arima-bonas-temp-residuals.dat" )
#dput( actual.wind-pred.wind[,1], "arima-bonas-wind-residuals.dat" )


