#Requirements
install.packages(c("dplyr", "tseries", "forecast", "zoo"))
library(dplyr)
library(zoo)
library(forecast)
library(tseries)

# Install and load the library for advanced ADF tests
if (!require("fUnitRoots")) install.packages("fUnitRoots")
library(fUnitRoots)

# Read CSV file manually (without readr or other libraries)
data_raw <- read.csv(
  "valeurs_mensuelles.csv",
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

### === ADF + Ljung-Box Combined Stationarity Check === ###
# Function to test white noise on ADF residuals
Qtests <- function(series, k = 24, fitdf = 0) {
  pvals <- sapply(1:k, function(l) {
    if (l <= fitdf) return(NA)
    Box.test(series, lag = l, type = "Ljung-Box", fitdf = fitdf)$p.value
  })
  return(data.frame(lag = 1:k, pval = pvals))
}

# Function to find ADF lag length with valid residuals (white noise)
adfTest_valid <- function(series, kmax = 24, type = "nc") {
  k <- 0
  repeat {
    cat(paste0("ADF with ", k, " lags: residuals OK? "))
    test <- suppressWarnings(adfTest(series, lags = k, type = type))
    fitdf <- length(test@test$lm$coefficients)
    pvals <- Qtests(test@test$lm$residuals, 24, fitdf)$pval
    if (all(is.na(pvals) | pvals > 0.05)) {
      cat("✔️ Residuals are white noise.\n")
      return(test)
    } else {
      cat("❌ Autocorrelation detected.\n")
      k <- k + 1
      if (k > kmax) stop("No valid ADF test found within lag limit.")
    }
  }
}

# Apply to each transformation of the series
adf_original <- adfTest_valid(ts_data, type = "ct")
adf_log <- adfTest_valid(log_ts, type = "ct")
adf_diff <- adfTest_valid(diff_ts, type = "nc")
adf_logdiff <- adfTest_valid(logdiff_ts, type = "nc")

# Summarize results
cat("\n--- ADF Test Results (test statistic) ---\n")
print(c(Original = adf_original@test$p.value,
        Log = adf_log@test$p.value,
        Diff = adf_diff@test$p.value,
        LogDiff = adf_logdiff@test$p.value))

cat("\n--- Ljung-Box on ADF Residuals (LogDiff series) ---\n")
print(Qtests(adf_logdiff@test$lm$residuals, 24))

#PART II: 
# Plot ACF and PACF to check autocorrelation
par(mfrow = c(1, 2))

acf(logdiff_ts, main = "ACF of Log-Differenced Series")
pacf(logdiff_ts, main = "PACF of Log-Differenced Series")

# We will test multiple ARMA models based on the ACF and PACF analysis. 

# Testing MA(2)
ma2_model <- arima(logdiff_ts, order = c(0, 0, 2))
summary(ma2_model)

# Testing ARMA(1,2)
arma12_model <- arima(logdiff_ts, order = c(1, 0, 2))
summary(arma12_model)

# Testing ARMA(2,2)
arma22_model <- arima(logdiff_ts, order = c(2, 0, 2))
summary(arma22_model)

# Testing ARMA(4,2)
arma42_model <- arima(logdiff_ts, order = c(4, 0, 2))
summary(arma42_model)

# Checking residuals for all models
# We will use the Ljung-Box test to assess if the residuals resemble white noise.
checkresiduals(ma2_model)
checkresiduals(arma12_model)
checkresiduals(arma22_model)
checkresiduals(arma42_model)

#Forecast
# load the necessary libraries
library(forecast)

# Forecast the next two periods (t+1, t+2)
forecast_values <- forecast(arma12_model, h = 2)

# Plot the forecast with x-axis starting from 2020
par(mfrow = c(1, 1))
plot(forecast_values, xlim = c(2020, max(time(forecast_values$mean))), main = "Forecast for t+1 and t+2", ylab = "Value", xlab = "Year")




# Plotting the original series with confidence forecast

forecast_values <- forecast(arma12_model, h = 2)

# Reverse the log-difference transformations
reverse_log_diff <- function(forecast_obj, last_value) {
  # Point forecasts
  point_fc <- numeric(length(forecast_obj$mean))
  point_fc[1] <- last_value * exp(forecast_obj$mean[1])
  point_fc[2] <- point_fc[1] * exp(forecast_obj$mean[2])
  
  # Lower CI (using the first confidence level)
  lower_fc <- numeric(length(forecast_obj$lower[,1]))
  lower_fc[1] <- last_value * exp(forecast_obj$lower[1,1])
  lower_fc[2] <- lower_fc[1] * exp(forecast_obj$lower[2,1])
  
  # Upper CI
  upper_fc <- numeric(length(forecast_obj$upper[,1]))
  upper_fc[1] <- last_value * exp(forecast_obj$upper[1,1])
  upper_fc[2] <- upper_fc[1] * exp(forecast_obj$upper[2,1])
  
  return(list(mean = point_fc, lower = lower_fc, upper = upper_fc))
}

last_original <- tail(ts_data, 1)  # Last value of original series
print(last_original)

reversed <- reverse_log_diff(forecast_values, last_original)
print(reversed)

# Create time series objects for the forecast period
fc_dates <- seq(end(ts_data), by = deltat(ts_data), length.out = 3)[-1]
fc_series <- zoo(reversed$mean, fc_dates)
lower_series <- zoo(reversed$lower, fc_dates)
upper_series <- zoo(reversed$upper, fc_dates)

# Combine with original series
combined <- merge(ts_data = ts_data, 
                 forecast = fc_series, 
                 lower = lower_series, 
                 upper = upper_series)

# Plotting everything
plot(combined$ts_data, 
     main = "Original Series with Forecasts", 
     xlab = "Time", ylab = "Value",
     ylim = range(combined, na.rm = TRUE),
     xlim = c(as.yearmon("2020-01"), end(ts_data)),
     lwd = 2)

lines(combined$forecast, col = "red", lwd = 2, type = "o", pch = 19)

# Add confidence interval
lines(combined$lower, col = "blue", lty = 2)
lines(combined$upper, col = "blue", lty = 2)

# Shaded confidence area
x_points <- c(
  as.numeric(fc_dates[1]),  # lower_1 x
  as.numeric(fc_dates[2]),  # lower_2 x
  as.numeric(fc_dates[2]),  # upper_2 x
  as.numeric(fc_dates[1])   # upper_1 x
)

y_points <- c(
  reversed$lower[1],  # lower_1 y
  reversed$lower[2],  # lower_2 y
  reversed$upper[2],  # upper_2 y
  reversed$upper[1]   # upper_1 y
)

polygon(x_points, y_points, col = rgb(0, 0, 1, 0.2), border = NA)

#Drawing the line to link the original series and the forecast
last_point_time <- time(tail(ts_data, 1))
first_forecast_point_time <- time(head(fc_series, 1))
lines(c(last_point_time, first_forecast_point_time),
      c(coredata(tail(ts_data, 1)), coredata(head(fc_series, 1))),
      col = "red", lwd = 2)

# Legend
legend("topleft", 
       legend = c("Original", "Forecast", "95% CI"), 
       col = c("black", "red", "blue"), 
       lty = c(1, 1, 2), 
       lwd = c(2, 2, 1),
       pch = c(NA, 19, NA))
