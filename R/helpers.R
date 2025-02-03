
#
# auxiliary functions
#
# (C) 2025 A. Steland, Institute of Statistics and AI Center, RWTH Aachen University
# published under the Academic Research Software Licence (ARSL) 1.0
#


MSE = function(x,y) as.numeric( crossprod(x-y) )
RMSE = function(x,y) as.numeric( sqrt(crossprod(x-y)) )

add_nonlinearities = function(dat) {
  
  dat$wind.res=Lag(dat$windspeedmph,-1)
  dat$tmpf.res=Lag(dat$tempf,-1)
  
  dat$solarradiation = c( 0, diff(dat$solarradiation))
  
  # # for best wind model exclude humidity and add lagged wind
  # dat = dat[,c(-6)]
  dat$wind.l1=Lag(dat$windspeedmph,1)
  dat$wind.l2=Lag(dat$windspeedmph,2)
  dat$wind.l3=Lag(dat$windspeedmph,3)
  dat$windroot10 = (dat$windspeedmph)^0.1
  
  dat$temp.l1=Lag(dat$tempf,1)
  dat$temp.l2=Lag(dat$tempf,2)
  dat$temp.l3=Lag(dat$tempf,3)
  dat$temproot10 = (dat$tempf)^0.1
  
  # nonlinear terms: quadratic and logs, increments and squared increments
  
  dat$wind2 = dat$windspeedmph^2
  dat$wind.log = log(dat$windspeedmph+1)
  dat$wind_delta = c(0, diff( dat$windspeedmph )) # /(dat$windspeedmph+1)
  dat$wind_delta_sq = c(0, diff( dat$windspeedmph ))^2
  
  
  dat$tempf2 = dat$tempf^2
  dat$tempf3 = dat$tempf^3
  dat$temp.log = log(dat$tempf+1)
  dat$temp_delta = c(0, diff( dat$tempf )) # /(dat$tempf)
  dat$temp_delta_sq = c(0, diff( dat$tempf ))^2
  
  # threshold effects
  
  dat$tempf.M = as.numeric( dat$tempf > quantile(dat$tempf, 0.5 ) )
  dat$tempf.Q3 = as.numeric( dat$tempf > quantile(dat$tempf, 0.75 ) )
  dat$windt.M = as.numeric( dat$tempf > quantile(dat$tempf, 0.5 ) )
  dat$windt.Q3 = as.numeric( dat$tempf > quantile(dat$tempf, 0.75 ) )
  dat$windt.Tail = as.numeric( dat$windspeedmph > quantile(dat$windspeedmph, 0.9 ) )
  
  dat
}

generate_table = function( actual.temp, actual.wind, pred.temp, pred.wind ) {
  
  MSE.wind = MSE( actual.wind, pred.wind[,1] )/nrow(pred.wind)
  MSE.wind
  
  MSE.temp = MSE( actual.temp, pred.temp[,1] )/nrow(pred.temp)
  MSE.temp
  
  # coverage probs
  CP.temp = mean((actual.temp>pred.temp[,2])*(actual.temp<pred.temp[,3]))
  CP.wind = mean((actual.wind>pred.wind[,2])*(actual.wind<pred.wind[,3]))
  
#  IS.x  = mean( wis( observed=actual, predicted=cbind(pred.x[,2], pred.x[,3]), quantile_level=c(0.025, 0.975) ) )
  IS.temp  = mean( wis( observed=actual.temp, predicted=cbind(pred.temp[,2], pred.temp[,3]), quantile_level=c(0.025, 0.975) ) )
  IS.wind  = mean( wis( observed=actual.wind, predicted=cbind(pred.wind[,2], pred.wind[,3]), quantile_level=c(0.025, 0.975) ) )
#  IS.temp  = mean( interval_score( actual.temp, pred.temp[,2], pred.temp[,3], 95  ) )
#  IS.wind = mean( interval_score( actual.wind, pred.wind[,2], pred.wind[,3], 95  ) )
  
  CRPS.temp = mean(crps_sample(actual.temp,as.matrix(pred.temp[,1])))
  CRPS.wind = mean(crps_sample(actual.wind,as.matrix(pred.wind[,1])))
  
  tab = rbind( c(MSE.temp, CRPS.temp, CP.temp, IS.temp ), 
               c(MSE.wind, CRPS.wind, CP.wind, IS.wind ) )
  colnames(tab) = c("MSE", "CRPS", "95% Coverage", "Interval score")
  rownames(tab) = c("Temp", "Wind" )
  tab
}

# Long-Run-Variance Estimator using Bartlett weights

lrvest = function( x, l ) {
  n = length(x)
  x = x - mean(x)
  l = floor( 12 * exp( 0.25 * log(n/100) ) )
  s = (n-1)/n * var(x)
  for( h in (1:l) ) {
    s = s + 2 * (1-h/(l+1)) * crossprod( x[1:(n-h)], x[(1+h):n] ) / n
  }
  s
}

#
# CI for E(MSE) based on given residuals (i.e. MSE = mean(residuals^2)). 
# Asymptotic variance estimated by Newey-West long-run-variance estimator
# 

CI_MSE = function( residuals, alpha = 0.05 ) {
  n = length(residuals) 
  MSE.hat = mean( residuals^2 )
  lrv.est = lrvest( residuals^2 )
  l = MSE.hat + sqrt(lrv.est) * qnorm(alpha/2) / sqrt(n)
  u = MSE.hat + sqrt(lrv.est) * qnorm(1-alpha/2) / sqrt(n)
  cat( "lower CI bound for MSE: ", l, "\n" )
  cat( "upper CI bound for MSE: ", u, "\n" )
  list( MSE = MSE.hat, lower = l, upper = u )
}

#
# Compare two MSE's using the fact that MSE1 - MSE2 = (1/n) sum_i (e_i^2 - f_i^2)
# where e_i are the residuals of model 1 and f_i the residuals of model 2,
# assuming they are calculated from the same sample of size n
# Asymptotic variance estimated by Newey-West long-run-variance estimator
# 


Test_Models_MSE = function( residuals1, residuals2, alpha = 0.05 ) {
  n = length(residuals1) 
  delta.hat = mean( residuals1^2 - residuals2^2 )
  lrv.est = lrvest( residuals1^2 - residuals2^2 )
  l = delta.hat + sqrt(lrv.est) * qnorm(alpha/2) / sqrt(n)
  u = delta.hat + sqrt(lrv.est) * qnorm(1-alpha/2) / sqrt(n)
  p = 2 * ( 1-pnorm( sqrt(n) * abs(delta.hat)/ sqrt(lrv.est) ) )
  cat("MSE1 = ", mean(residuals1^2), "\n")
  cat("MSE2 = ", mean(residuals2^2), "\n")
  cat("Delta = ", mean(residuals1^2-residuals2^2), "\n")
  cat( "lower CI bound for difference in expected MSE1 - MSE2: ", l, "\n" )
  cat( "upper CI bound for difference in expected MSE1 - MSE2: ", u, "\n" )
  list( delta = delta.hat, lower = l, upper = u, p = p )
}

