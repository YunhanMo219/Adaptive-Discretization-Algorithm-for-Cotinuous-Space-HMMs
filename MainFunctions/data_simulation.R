# There are two functions in data_simulation.R. "sim_count_hmm" and "plot_poisson_series"

library(data.table)
library(ggplot2)

# define data simulation function in 1d setting.
# hidden state: AR(1), noise~normal;
# observed Y ~ poisson(w*beta*e^X_t)
sim_count_hmm <- function(TT, phi, sigma, beta, w){
  # set.seed(123)
  
  DT <- data.table(
    time = 1:TT,
    X = numeric(TT),
    lambda = numeric(TT),
    Y = numeric(TT)
  )
  
  DT$X[1] <- rnorm(1, mean = 0, sd = sqrt(sigma^2 / (1-phi^2)))
  DT$lambda[1] <- w*beta*exp(DT$X[1])
  DT$Y[1] <- rpois(1, lambda = DT$lambda[1])
  
  for (t in 2:TT) {
    DT$X[t] <- rnorm(1, mean = phi*DT$X[t-1], sd = sigma)
    DT$lambda[t] <- w*beta*exp(DT$X[t])
    DT$Y[t] <- rpois(1, lambda = DT$lambda[t])
  }
  return(DT)
}

# plot for observed data Y
plot_poisson_series <- function(data, phi, sigma, beta) {
  
  plot_title <- sprintf("Poisson Time Series (phi = %.2f, sigma = %.2f, beta = %.2f)",
                        phi, sigma, beta)
  
  p <- ggplot(data, aes(x = time, y = Y)) +
    geom_line(color = "black") +
    geom_point(color = "black", size = 0.9) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Y"
    ) +
    theme_minimal()
  
  print(p)
}


# Function to plot true vs estimated hidden state using ggplot2
plot_hidden_state_series <- function(true_state, estimated_state, phi, sigma, beta) {
  
    df <- data.frame(
    Time = seq_along(true_state),
    TrueState = true_state,
    EstimatedState = estimated_state
  )
  
  # Reshape to long format for ggplot
  df_long <- tidyr::pivot_longer(
    df,
    cols = c("TrueState", "EstimatedState"),
    names_to = "StateType",
    values_to = "Value"
  )
  
  MSE <- sum((true_state - estimated_state)^2) / length(true_state)
  
  plot_title <- sprintf(
    "Hidden State: phi = %.2f, sigma = %.2f, beta = %.2f, MSE = %.3f",
    phi, sigma, beta, MSE
  )
  
  
  p <- ggplot(df_long, aes(x = Time, y = Value, color = StateType)) +
    geom_line(size = 0.6) +
    scale_color_manual(
      values = c("EstimatedState" = "black", "TrueState" = "red"),
      labels = c("Estimated Hidden State", "Real Hidden State"),
      name = NULL
    ) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Hidden State"
    ) +
    theme_minimal() 
  
    print(p)
}


