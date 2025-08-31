# read data
eq6 <- read.csv("earthquake6_yearly_count.csv")

# source code
base_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"
source(file.path(base_dir, "auto_functions.R"))
source(file.path(base_dir, "useful_functions.R")) # adapted from Zimmerman et al (2024)
source(file.path(base_dir, "HMM_Functions.R")) # from Zimmerman et al (2024)
base_dir2 <- "E:/study in the UK/research project/My Code/Simulation3"
source(file.path(base_dir2, "fixed discretization function.R"))
base_dir3 <- "E:/study in the UK/research project/My Code/Earthquake"
source(file.path(base_dir3, "hmm_prediction.R"))

# copy the estimation from adaptive algorithm, as Fridman and Harris
# suggest using −cmin = cmax = 3σx, where σx denotes the stationary 
# standard deviation of Xt
# eq_beta_hat <- 104.3305
# eq_phi_hat <- 0.9754318
# eq_sigma_hat <- 0.1768331

eq_beta_hat <- 104.2548
eq_phi_hat <- 0.9754318
eq_sigma_hat <- 0.1768331

sigma_x <- sqrt(eq_sigma_hat^2 / (1 - eq_phi_hat^2))
b <- 3 * sigma_x
m <- 100
w <- 1
Y <- eq6$count

fixed_results <- local_decoding_fixed(Y = Y, w = w, b = b, m = m)

# local decoding mse
fixed_mse <- mean((fixed_results$lambda_hat - eq6$count)^2)
fixed_mse

# probability of true count fall into 95% confidence interval
fix_decoding_converged <- mean(eq6$count >= fixed_results$lambda_lower_eti &  
                                eq6$count <= fixed_results$lambda_upper_eti)
fix_decoding_converged


# window prediction
window_size <- 100          # Use 100 observations each time
forecast_horizon <- 1       # Forecast 1 step ahead each time
x_list <- 10:260             # Possible forecasted counts

# Initialize list to store forecast distributions
forecast_list <- list()

# Loop index: from 1 to (total length - window size)
for (i in 1:(length(Y) - window_size)) {
  
  # Extract current rolling window
  Y_window <- Y[i:(i + window_size - 1)]
  
  # Estimate model and do local decoding
  results <- local_decoding_fixed(Y = Y_window, w = w, b = b, m = m)
  
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
  file = "earthquake6_forecast_all_fixed.csv",
  row.names = TRUE  # Keep row names as x_list
)


# get prediction
pred_prob <- read.csv("earthquake6_forecast_all_fixed.csv")
y_predict_fixed <- sapply(pred_prob[,-1], function(prob_y){
  id <- which.max(prob_y)
  pred_prob[,1][id]
})

ci_95 <- sapply(pred_prob[,-1], function(prob_y){
  cdf <- cumsum(prob_y)
  lower_id <- which(cdf >= 0.025)[1]
  upper_id <- which(cdf >= 0.975)[1]
  c(pred_prob[lower_id, 1], pred_prob[upper_id, 1])
})
ci_95 <- t(ci_95)
colnames(ci_95) <- c("Lower95", "Upper95")

# Compute coverage of 95% equal-tailed intervals
true_values <- Y[101:125]
coverage_fixed <- mean(true_values >= ci_95[, "Lower95"] & 
                         true_values <= ci_95[, "Upper95"])
coverage_fixed

# Compute MSE between observed and predicted counts
mse <- mean((true_values - y_predict_fixed)^2)
mse

# plots for check
num_plots <- ncol(forecast_matrix)
x_vals <- as.numeric(rownames(forecast_matrix))
years <- 2000:2024
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

# pseudo residual test
pseudo_residual <- numeric()
for (i in 1:length(true_values)) {
  P <- sum(pred_prob[,i+1][1:true_values[i]])
  pseudo_residual[i] <- qnorm(P)
}

plot(1:25, pseudo_residual)
shapiro.test(pseudo_residual)
pseudo_residual


