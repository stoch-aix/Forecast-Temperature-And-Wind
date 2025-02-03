
#
# Template from Gemini 2.0 flash
#


# Install packages if not already installed
if (!require(quantregForest)) install.packages("quantregForest")
if (!require(caret)) install.packages("caret")
if (!require(dplyr)) install.packages("dplyr")

# Load libraries
library(quantregForest)
library(caret)
library(dplyr)

# Your data frame
# data <- read.csv("your_data.csv")  # Uncomment and modify to load your data

# Example structure of data
 data <- data.frame(
   Date = as.Date('2021-01-01') + 0:364,
   Temperature = runif(365, min = -10, max = 35),
   WindSpeed = runif(365, min = 0, max = 15),
   SolarIrradiance = runif(365, min = 0, max = 1000),
   Humidity = runif(365, min = 10, max = 100)
 )

# Create lagged variables (lag of 1 day)
data <- data %>%
  arrange(Date) %>%
  mutate(
    Temp_Lag1 = lag(Temperature, 1),
    WindSpeed_Lag1 = lag(WindSpeed, 1),
    SolarIrradiance_Lag1 = lag(SolarIrradiance, 1),
    Humidity_Lag1 = lag(Humidity, 1)
  )

# Remove rows with NA values resulting from lagging
data <- na.omit(data)

# Set seed for reproducibility
set.seed(123)

# Split data: 80% training, 20% testing
trainIndex <- createDataPartition(data$Temperature, p = 0.8, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

# Define the model formula
model_formula <- Temperature ~ Temp_Lag1 + WindSpeed_Lag1 + SolarIrradiance_Lag1 + Humidity_Lag1

# Train the Quantile Regression Forest model
qrf_model <- quantregForest(
  x = trainData %>% select(Temp_Lag1, WindSpeed_Lag1, SolarIrradiance_Lag1, Humidity_Lag1),
  y = trainData$Temperature,
  ntree = 1000  # Number of trees; increase for better performance
)


# Prepare test predictors
testPredictors <- testData %>% select(Temp_Lag1, WindSpeed_Lag1, SolarIrradiance_Lag1, Humidity_Lag1)

# Predict the mean temperature
mean_prediction <- predict(qrf_model, newdata = testPredictors)

# Predict the 2.5% and 97.5% quantiles for prediction intervals
quantiles <- predict(qrf_model, newdata = testPredictors, what = c(0.025, 0.975))

# Combine predictions and intervals
results <- testData %>%
  select(Date, Actual_Temperature = Temperature) %>%
  mutate(
    Predicted_Temperature = mean_prediction,
    Lower_Bound = quantiles[, 1],
    Upper_Bound = quantiles[, 2]
  )

# View the first few results
head(results)

# Calculate prediction accuracy metrics
library(Metrics)
mae <- mae(results$Actual_Temperature, results$Predicted_Temperature)
rmse <- rmse(results$Actual_Temperature, results$Predicted_Temperature)

cat("Mean Absolute Error:", round(mae, 2), "\n")
cat("Root Mean Squared Error:", round(rmse, 2), "\n")

plot( results$Date, results$Actual_Temperature, type="l", ylim=c(-15,45) )
lines( results$Date, results$Predicted_Temperature[,2], col="blue")
lines( results$Date, results$Lower_Bound, col="blue", lty=2)
lines( results$Date, results$Upper_Bound, col="blue", lty=2)

# Plot actual vs predicted temperatures with prediction intervals
library(ggplot2)

p <- ggplot(results, aes(x = Date)) 

p <- p + geom_line(aes(y = Actual_Temperature), size = 1, linetype = "dashed") 

p <- p + geom_line(aes(y = Predicted_Temperature), size = 1)


  geom_ribbon(aes(ymin = Lower_Bound, ymax = Upper_Bound), alpha = 0.2) +
  labs(
    title = "Actual vs Predicted Temperature",
    y = "Temperature",
    x = "Date"
  )
#+  theme_minimal()
