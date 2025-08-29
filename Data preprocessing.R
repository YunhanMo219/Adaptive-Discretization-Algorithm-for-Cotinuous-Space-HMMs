library(tidyverse)
library(lubridate)

code_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"

source(file.path(code_dir,"data_simulation.R"))
source(file.path(code_dir,"HMM_Functions.R")) # from Zimmerman et al (2024)
source(file.path(code_dir,"useful_functions.R")) # adapted from Zimmerman et al (2024)
source(file.path(code_dir,"auto_functions.R"))

data_dir <- "E:/study in the UK/research project/My Code/Earthquake"
data <- read.csv(file.path(data_dir, "earthquake6.csv"))
data$time <- as.POSIXct(data$time, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
data$year <- format(data$time, "%Y")

year_counts <- table(data$year)
print(year_counts)

year_counts_df <- as.data.frame(year_counts)
names(year_counts_df) <- c("year", "count")


year_counts_df$year <- as.numeric(as.character(year_counts_df$year))
write.csv(year_counts_df, "earthquake6_yearly_count.csv", row.names = FALSE)

# 画折线图
plot(year_counts_df$year, year_counts_df$count, 
     type = "o", 
     pch = 20,
     cex = 0.7,
     xlab = "Year", 
     ylab = "Count", 
     main = "Number of Earthquakes in the World")


Y <- as.matrix(year_counts_df$count)
w <- 1
results <- local_decoding_auto(Y,w)

saveRDS(results, file = "earthquake6_results.rds")
result_list <- readRDS("earthquake6_results.rds")

plot(results$Xt_hat,
     type = 'l',
     xlab = "Year", 
     ylab = "Hidden State", 
     main = 'Estimsted Hidden State')

# Parametric Bootstrap Prediction Interval
lambda_hat <- result_list$para_est$beta_hat * w * exp(result_list$Xt_hat)
n_B <- 100
TT <- length(lambda_hat)

sim_vector <- rpois(n = TT * n_B, lambda = rep(lambda_hat, each = n_B))

sim_matrix <- matrix(sim_vector, nrow = TT, ncol = n_B, byrow = TRUE)

lower <- apply(sim_matrix, 1, quantile, probs = 0.025)
upper <- apply(sim_matrix, 1, quantile, probs = 0.975)

inside <- as.numeric(Y) >= lower & as.numeric(Y) <= upper


plot_df <- tibble(
  time = year_counts_df$year,
  obs = as.numeric(Y),
  lower = lower,
  upper = upper,
  inside = inside  # TRUE / FALSE
)

ggplot(plot_df, aes(x = time)) +
  geom_linerange(aes(ymin = lower, ymax = upper), color = "black") +
  geom_point(aes(y = obs, color = inside), size = 2) +
  labs(
    title = "95% Parametric Bootstrap Prediction Interval & Observed Counts",
    x = "Year",
    y = "Count",
    color = "Inside 95% Interval"
  ) +
  theme_minimal()

## Plot for lambda hat
# Create a combined data frame for plotting
plot_df <- data.frame(
  year = year_counts_df$year,
  count = year_counts_df$count,
  lambda_hat = results$lambda_hat,
  lower = results$lambda_lower_eti,
  upper = results$lambda_upper_eti
)


ggplot(plot_df, aes(x = year)) +
  # 95% prediction band
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% ETI"), alpha = 0.4) +
  # Observed: line + points (solid circle)
  geom_line(aes(y = count, color = "Observed"), linewidth = 0.7) +
  geom_point(aes(y = count, color = "Observed", shape = "Observed"), size = 2) +
  # Estimated: line + points (solid triangle)
  geom_line(aes(y = lambda_hat, color = "Estimated"), linewidth = 0.7) +
  geom_point(aes(y = lambda_hat, color = "Estimated", shape = "Estimated"), size = 2) +
  scale_color_manual(
    name = NULL,
    values = c("Observed" = "black", "Estimated" = "red")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Observed" = 16, "Estimated" = 17)  # 16 = circle, 17 = triangle
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("95% ETI" = "orange")
  ) +
  labs(
    y = "Count",
    title = "Estimated Number of Earthquakes (>6.0) from Local Decoding"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Mean squared error between estimated and true counts
mse <- mean((plot_df$lambda_hat - plot_df$count)^2)
mse
# Coverage: probability that lambda_hat falls inside 95% equal-tailed interval
coverage <- mean(plot_df$count >= plot_df$lower & plot_df$count <= plot_df$upper)
coverage



