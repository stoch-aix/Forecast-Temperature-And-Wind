#
# Nonlinear robust ARIMA prediction of temperature and wind
# 
# Steland, A. (2024). Discussion on Assessing predictability of environmental time series with statistical and machine learning models  
# 
# (C) A. Steland, Institute of Statistics and AI Center, RWTH Aachen University
# published under the Academic Research Software Licence (ARSL) 1.0
#


library(forecast)
library(Hmisc)
library(dplyr)
library(fastDummies)

library(scoringutils)

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

# discard first rows due to NAs (lagged values)

train = train[4:nrow(train),]

# MSE = 0.467
regressors.wind = c("tempf", "temp_delta", "windt.Tail", "wind2", "wind.log", "wind_delta", "wind_delta_sq" )

# Alternative model using only tempf, MSE = 0.484
#regressors.wind = c( "tempf", "windt.Q3", "windroot10", "wind.log", "wind_delta", "wind_delta_sq" )

xreg.wind.train=data.matrix(train[,c(regressors.wind)])

regressors.temp = c("humidity", "solarradiation", "tempf2", "temp.log", 
                    "temp_delta", "temp_delta_sq")
xreg.temp.train=data.matrix(train[,regressors.temp])

# 
xreg.wind.test=test[,regressors.wind]
xreg.temp.test=test[,regressors.temp]

# modified forecast procedure
# using robust White estimator of the variance-covariance matrix of the OLS estimator

myforecast = function( model, xreg, h = 1, bootstrap=F, xreg.train = NULL ) {
  res = forecast( object=model, xreg=xreg, h=h, bootstrap=bootstrap ) 
  # change notation...
  xnew = xreg
  xreg = as.matrix( xreg.train )
  Xm = t( xreg ) %*% xreg / nrow(xreg)
  Dm = t( xreg ) %*% diag( model$residuals^2 ) %*% xreg / nrow(xreg)
  Xm.inv = solve(Xm)
  SWhite = Xm.inv %*% solve(Dm) %*% Xm.inv
  q.8 = qnorm(0.9)
  q.95 = qnorm(0.975)
  for ( i in (1:nrow(xnew) ) ) {
    v = sqrt( t(xnew[i,]) %*% SWhite %*% xnew[i,] / nrow(xreg) + var(model$residuals) ) 
    res$lower[i,1] = res$mean[i] - q.8 * v
    res$lower[i,2] = res$mean[i] - q.95 * v
    res$upper[i,1] = res$mean[i] + q.8 * v
    res$upper[i,2] = res$mean[i] + q.95 * v
  }
  res
}


# ARIMA

wind.model=auto.arima(train$wind.res,xreg = xreg.wind.train, seasonal = FALSE)
temp.model=auto.arima(train$tmpf.res,xreg = xreg.temp.train, seasonal = FALSE)

pred.wind=matrix(nrow=(nrow(test)-1),ncol=3)
pred.temp=matrix(nrow=(nrow(test)-1),ncol=3)

wind.model2=Arima(train$wind.res,xreg=as.matrix(xreg.wind.train),model=wind.model)

#wind.1=data.frame(forecast(wind.model2,xreg=as.matrix(xreg.wind.test[1,]),h=1))
wind.1=data.frame(myforecast(wind.model2,xreg=as.matrix(xreg.wind.test[1,]),h=1,xreg.train=xreg.wind.train))
pred.wind[1,]=c(wind.1$Point.Forecast,wind.1$Lo.95,wind.1$Hi.95)

temp.model2=Arima(train$tmpf.res,xreg=as.matrix(xreg.temp.train),model=temp.model)

#temp.1=data.frame(forecast(temp.model2,xreg=as.matrix(xreg.temp.test[1,]),h=1))
temp.1=data.frame(myforecast(temp.model2,xreg=as.matrix(xreg.temp.test[1,]),h=1,xreg.train=xreg.temp.train))
pred.temp[1,]=c(temp.1$Point.Forecast,temp.1$Lo.95,temp.1$Hi.95)

### Forecast
for (i in 2:nrow(pred.wind)){
  print(i)
  wind.model2=Arima(c(train$wind.res,test$wind.res[1:(i-1)]),xreg=as.matrix(rbind(xreg.wind.train,xreg.wind.test[1:(i-1),])),model=wind.model)
  #wind.1=data.frame(forecast(wind.model2,xreg=as.matrix(xreg.wind.test[i,]),h=1,bootstrap=F))
  wind.1=data.frame(myforecast(wind.model2,xreg=as.matrix(xreg.wind.test[i,]),h=1,xreg.train=as.matrix(rbind(xreg.wind.train,xreg.wind.test[1:(i-1),]))))
  pred.wind[i,]=c(wind.1$Point.Forecast,wind.1$Lo.95,wind.1$Hi.95)
  
  temp.model2=Arima(c(train$tmpf.res,test$tmpf.res[1:(i-1)]),xreg=as.matrix(rbind(xreg.temp.train,xreg.temp.test[1:(i-1),])),model=temp.model)
  #temp.1=data.frame(forecast(temp.model2,xreg=as.matrix(xreg.temp.test[i,]),h=1,bootstrap=F))
  temp.1=data.frame(myforecast(temp.model2,xreg=as.matrix(xreg.temp.test[i,]),h=1,xreg.train=as.matrix(rbind(xreg.temp.train,xreg.temp.test[1:(i-1),]))))
  pred.temp[i,]=c(temp.1$Point.Forecast,temp.1$Lo.95,temp.1$Hi.95)
  
}


### Actual values
actual.temp=test$tmpf.res[1:(nrow(pred.temp)-1)]
actual.wind=test$wind.res[1:(nrow(pred.wind)-1)]

### Replace predicted negative values with 0
pred.wind[pred.wind<0]<-0

tab = generate_table( actual.temp, actual.wind, 
                      pred.temp[1:(nrow(pred.temp)-1),], pred.wind[1:(nrow(pred.wind)-1),] )
tab

sub = 1:50
plot((actual.temp[sub]-32)/1.8,type="l", main="Temperature", ylim=c(5,30))
lines((pred.temp[sub,1]-32)/1.8, col="blue" )

#sub = 63:(63+28)
plot( actual.wind[sub],type="l", main="Wind" )
lines( pred.wind[sub,1], col="blue" )

#
# CI
# 

CI_MSE( actual.temp-pred.temp[,1] )

dput( actual.temp-pred.temp[,1], "arima-temp-residuals.dat" )
dput( actual.wind-pred.wind[,1], "arima-wind-residuals.dat" )

# 
# Test stationarity
# 

# KPSS test der Residuen akzeptiert H0: Stationaritaet
kpss.test( temp.model$residuals, null=c("Level") )
kpss.test( wind.model$residuals, null=c("Level") )




