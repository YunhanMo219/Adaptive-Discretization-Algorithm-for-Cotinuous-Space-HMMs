library(ggplot2)
library(forecast)

# read results and data
pred_prob <- read.csv("earthquake6_forecast_all.csv")
year_counts_df <- read.csv("earthquake6_yearly_count.csv")

Y <- as.matrix(year_counts_df$count)
true_values <- Y[101:125]

# compute pseudo-residual
pseudo_residual <- numeric()
for (i in 1:length(true_values)) {
  P <- sum(pred_prob[,i+1][1:true_values[i]])
  pseudo_residual[i] <- qnorm(P)
}

plot(1:25, pseudo_residual)
shapiro.test(pseudo_residual)
pseudo_residual

# plot acf of pseudo residual
acf_obj <- acf(pseudo_residual, plot = FALSE)

df_acf <- data.frame(
  Lag = acf_obj$lag,
  ACF = acf_obj$acf
)

# acf plot
ggplot(df_acf, aes(x = Lag, y = ACF)) +
  geom_segment(aes(xend = Lag, yend = 0), color = "steelblue", linewidth = 0.6) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-1.96/sqrt(length(pseudo_residual)), 
                            1.96/sqrt(length(pseudo_residual))),
             linetype = "dashed", color = "red") +
  labs(title = "ACF of Pseudo Residuals") +
  theme_minimal(base_size = 12)

# plot predicted Y and 95% eti for 2000 to 2024
y_predict <- sapply(pred_prob[,-1], function(prob_y){
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

years <- 2000:2024

plot_df <- data.frame(
  Year     = years,
  Observed = true_values,
  Pred     = y_predict,
  Lower    = ci_95[, "Lower95"],
  Upper    = ci_95[, "Upper95"]
)


ggplot(plot_df, aes(x = Year)) +
  # 95% prediction band
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "95% ETI"), alpha = 0.18) +
  # Observed: line + points (solid circle)
  geom_line(aes(y = Observed, color = "Observed"), linewidth = 0.7) +
  geom_point(aes(y = Observed, color = "Observed", shape = "Observed"), size = 2) +
  # Predicted (MAP): line + points (star)
  geom_line(aes(y = Pred, color = "Predicted"), linewidth = 0.7) +
  geom_point(aes(y = Pred, color = "Predicted", shape = "Predicted"), size = 2) +
  scale_color_manual(
    name = NULL,
    values = c("Observed" = "black", "Predicted" = "red")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Observed" = 16, "Predicted" = 17) 
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("95% ETI" = "red")
  ) +
  labs(
    y = "Count",
    title = "Predicted Number of Earthquakes (>6.0) Using 100-Year Rolling Window"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Compute coverage of 95% equal-tailed intervals
inside <- (true_values >= ci_95[, "Lower95"]) & (true_values <= ci_95[, "Upper95"]) # indicator if true value is inside interval
coverage <- mean(inside) # empirical coverage rate
coverage

# Compute MSE between observed and predicted counts
mse <- mean((true_values - y_predict)^2)
mse

