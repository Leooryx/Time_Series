#Léo DONZIL and Léo LEROY

#Requirements
install.packages(c("dplyr", "tseries", "forecast", "zoo"))
library(dplyr)
library(zoo)
library(forecast)
library(tseries)

# Install and load the library for advanced ADF tests
if (!require("fUnitRoots")) install.packages("fUnitRoots")
library(fUnitRoots)

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

###########
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
    test <- adfTest(series, lags = k, type = type)
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
print(c(Original = adf_original@test$statistic,
        Log = adf_log@test$statistic,
        Diff = adf_diff@test$statistic,
        LogDiff = adf_logdiff@test$statistic))

cat("\n--- Ljung-Box on ADF Residuals (LogDiff series) ---\n")
print(Qtests(adf_logdiff@test$lm$residuals, 24, fitdf = length(adf_logdiff@test$lm$coefficients)))
###########



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