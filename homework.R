#Léo DONZIL and Léo LEROY

#Requirements
install.packages(c("dplyr", "tseries", "forecast", "zoo"))
library(dplyr)
library(zoo)
library(forecast)
library(tseries)

#PART I: The Data
# Read CSV file manually (without readr or other libraries)
data_raw <- read.csv(
  "/home/onyxia/Time_Series/valeurs_mensuelles.csv",
  sep = ";",         # Use semicolon as separator
  header = FALSE,    # Skip the header row
  skip = 4,          # Skip the first 4 rows (metadata)
  stringsAsFactors = FALSE
)


# Clean and transform the data
data <- data_raw %>% 
  select(Date = V1, Index = V2) %>%  # Rename columns to Date and Index
  mutate(
    Date = as.yearmon(Date),        # parse "YYYY-MM" to zoo-compatible format
    Index = as.numeric(Index)       # ensure it's numeric
  ) %>% 
  filter(!is.na(Date), !is.na(Index))  # clean rows with parsing errors

# Create zoo time series
ts_data <- zoo(data$Index, order.by = data$Date)

# Plot original series using base R
plot(ts_data, type = "l", col = "blue", lwd = 2, 
     main = "Original Time Series – Industrial Production Index (Pharma)", 
     ylab = "Index (Base 100 in 2021)", 
     xlab = "Year")


#Checking stationarity
adf.test(data$Index)  # Augmented Dickey-Fuller test
#the p-value being 0.1057, we do not reject that the time series is non-stationary

#Trying transformations to see improvements on the ADF test
# 1. Log transformation
log_ts <- log(ts_data)

# 2. First difference
diff_ts <- diff(ts_data)

# 3. Log + difference
logdiff_ts <- diff(log_ts)



par(mfrow = c(2, 2))

plot(ts_data, main = "Original Series", ylab = "Index")
plot(log_ts, main = "Log-transformed", ylab = "log(Index)")
plot(diff_ts, main = "Differenced", ylab = "Δ Index")
plot(logdiff_ts, main = "Log-Differenced", ylab = "Δ log(Index)")

#Checking stationarity for the different time series
check_stationarity <- function(x) {
  test <- adf.test(na.omit(x))
  return(c("ADF p-value" = test$p.value, "Variance" = var(na.omit(x))))
}

results <- rbind(
  Original = check_stationarity(ts_data),
  Log = check_stationarity(log_ts),
  Differenced = check_stationarity(diff_ts),
  LogDiff = check_stationarity(logdiff_ts)
)

print(results)
#            ADF p-value     Variance
#Original      0.7661196 5.693979e+02
#Log           0.6266858 1.466186e-01
#Differenced   0.0100000 1.922875e+01
#LogDiff       0.0100000 2.762987e-03
#Only Differenced and LogDiff are significantly not stationary.
#Since LogDiff has the smallest variance, we will work with this one



#PART II: 
# Plot ACF and PACF to check autocorrelation
par(mfrow = c(1, 2))

acf(logdiff_ts, main = "ACF of Log-Differenced Series")
pacf(logdiff_ts, main = "PACF of Log-Differenced Series")

# Auto ARIMA model selection
auto_model <- auto.arima(logdiff_ts)

# Display the chosen ARMA model
summary(auto_model)

# Check residuals for validity (should resemble white noise)
checkresiduals(auto_model)

# Part III: Prediction
# Forecast the next 12 months with 95% confidence level
forecast_horizon <- 12
pred <- forecast(auto_model, h = forecast_horizon, level = 95)

par(mfrow = c(1, 1))
# Plot the forecast with confidence intervals
plot(pred,
     main = "Forecast with 95% Confidence Interval",
     xlab = "Time",
     ylab = "Log-Differenced Index")

# Optional: Add grid for clarity
grid()