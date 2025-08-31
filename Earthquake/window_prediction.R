library(tidyverse)
library(lubridate)

code_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"

source(file.path(code_dir,"data_simulation.R"))
source(file.path(code_dir,"HMM_Functions.R")) # from Zimmerman et al (2024)
source(file.path(code_dir,"useful_functions.R")) # adapted from Zimmerman et al (2024)
source(file.path(code_dir,"auto_functions.R"))
base_dir3 <- "E:/study in the UK/research project/My Code/Earthquake"
source(file.path(base_dir3, "hmm_prediction.R"))


# --------------------------
# 1. Load data
# --------------------------
year_counts_df <- read.csv("earthquake6_yearly_count.csv")

# Convert to matrix or vector
Y <- as.matrix(year_counts_df$count)

# --------------------------
# 2. Parameters
# --------------------------
window_size <- 100          # Use 100 observations each time
forecast_horizon <- 1       # Forecast 1 step ahead each time
x_list <- 10:250             # Possible forecasted counts
w <- 1                      # Your parameter for local_decoding_auto()

# --------------------------
# 3. Loop through all possible windows and forecast
# --------------------------
# Initialize list to store forecast distributions
forecast_list <- list()

# Loop index: from 1 to (total length - window size)
for (i in 1:(length(Y) - window_size)) {
  
  # Extract current rolling window
  Y_window <- Y[i:(i + window_size - 1)]
  
  # Estimate model and do local decoding
  results <- local_decoding_auto(Y_window, w)
  
  # Extract estimated parameters
  para <- results$para_est
  ess_para <- list(
    cmin = results$cmin,
    cmax = results$cmax,
    m = results$m
  )
  
  # Forecast 1 step ahead
  forecast_1h <- poisson_hmm_forecast(
    para, ess_para, Y_window, w, x_list, forecast_horizon
  )
  
  # Save the forecasted probability distribution (only the first row)
  forecast_list[[i]] <- forecast_1h[1, ]
  
  # Optional: print progress
  cat("Done: window", i, "\n")
}

# --------------------------
# 4. Combine forecasts into a data frame
# --------------------------
# Convert list to matrix: rows = x_list, columns = each forecasted year
forecast_matrix <- do.call(cbind, forecast_list)

# Row names = possible count values
rownames(forecast_matrix) <- x_list

# Optional: column names = forecasted year index
colnames(forecast_matrix) <- paste0(
  "Year_", 2000:2024
)

# --------------------------
# 5. Save to CSV
# --------------------------
write.csv(
  forecast_matrix,
  file = "earthquake6_forecast_all.csv",
  row.names = TRUE  # Keep row names as x_list
)


# hist plot for each year
# --------------------------
# 1. Prepare data
# --------------------------
true_values <- Y[101:125]

# Years to annotate:
years <- 2000:2024

num_plots <- ncol(forecast_matrix)
x_vals <- as.numeric(rownames(forecast_matrix))

# --------------------------
# 2. Setup multi-panel plot
# --------------------------
par(mfrow = c(5, 5),  # 5 rows x 5 columns
    mar = c(3, 3, 2, 2))  # Margins: bottom, left, top, right

# --------------------------
# 3. Loop through each column and plot
# --------------------------
for (i in 1:num_plots) {
  
  probs <- forecast_matrix[, i]
  
  # X-axis = possible counts
  x <- x_vals
  
  # Y-axis = probability
  y <- probs
  
  # Plot empty plot
  plot(x, y, type = "n",
       xlab = "", ylab = "",
       main = "", axes = FALSE)
  
  # Draw blue sticks for the forecast distribution
  segments(x0 = x, y0 = 0, x1 = x, y1 = y, col = "lightblue")
  
  # Draw red stick for the true value
  true_val <- true_values[i]
  
  # Only draw if in range
  if (true_val >= min(x) && true_val <= max(x)) {
    true_idx <- which(x == true_val)
    segments(x0 = x[true_idx], y0 = 0, 
             x1 = x[true_idx], y1 = y[true_idx], 
             col = "red", lwd = 2)
  }
  
  # Add axes
  axis(1, at = seq(min(x), max(x), length.out = 5))
  axis(2, las = 1, cex.axis = 0.7)
  
  # Add year in top right corner
  usr <- par("usr")
  text(
    x = usr[2] * 0.95,
    y = usr[4] * 0.95,
    labels = years[i],
    adj = c(1, 1),
    cex = 1.0
  )
}

# --------------------------
# 4. Reset layout
# --------------------------
par(mfrow = c(1, 1))



