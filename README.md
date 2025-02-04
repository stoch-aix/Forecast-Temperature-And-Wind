## Prediction of Local Temperature and Wind Speed

This repository provides selected R source code used for the paper

​	*Steland, Ansgar (2025). Discussion on 'Assessing predictability of environmental time series with statistical and machine learning methods', Environmetrics*



**Files:**

`helpers.R:`

​	Auxiliary functions

`arima-bonas-reanalysis.R:`

​	Reproduces the results for the best ARIMA prediction approach and provides additional analyses.

`nlr-arima-public.R:`

​	Fits and analyses prediction models based on the nonlinear robust ARIMA approach discussed in the above discussion paper. These models outperform all machine learning methods investigated for this data set.

`Env-QuantRegForests-public.R`

​	Fits and analyses prediction based on quantile regression forests. 

`QuantRegForests-template.R`

​	Example code generated by the Gemini 2.0 flash LLN. Surprisingly correct. The corresponding prompt is provided in the above publication.

(C) 2025 Ansgar Steland, Institute of Statistics and AI Center, RWTH Aachen University
