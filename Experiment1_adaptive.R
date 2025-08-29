# Experiment 1 aims to compare the accuracy of parameter estimation between the 
# adaptive discretization algorithm and the fixed essential domain method. 
# In addition, we report the average length of the adaptive essential domain and
# the rule-of-thumb value suggested by Fridman and Harris.

# set of the simulation
beta <- 20
phi <- 0.9
# sigma <- 0.2 #Exp1
sigma <- 0.4
TT <- 1000

n_sim <- 200

# source functinos
simu_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"
source(file.path(simu_dir, "data_simulation.R"))

code_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"
source(file.path(code_dir,"data_simulation.R"))
source(file.path(code_dir,"HMM_Functions.R")) # from Zimmerman et al (2024)
source(file.path(code_dir,"useful_functions.R")) # adapted from Zimmerman et al (2024)
source(file.path(code_dir,"auto_functions.R"))

# Create output folder to store experiment 1 data
dir.create("Exp2_data", showWarnings = FALSE)

set.seed(123)
X_mat <- matrix(NA, nrow = TT, ncol = n_sim)
Y_mat <- matrix(NA, nrow = TT, ncol = n_sim)
lambda_mat <- matrix(NA, nrow = TT, ncol = n_sim)  # store lambda

for (i in 1:n_sim) {
  DT1 <- sim_count_hmm(TT, phi, sigma, beta, w=1)
  X_mat[, i] <- DT1$X
  Y_mat[, i] <- DT1$Y
  lambda_mat[, i] <- DT1$lambda
}

write.csv(X_mat, file.path("Exp2_data", "X_B_Sample.csv"), row.names = FALSE)
write.csv(Y_mat, file.path("Exp2_data", "Y_B_Sample.csv"), row.names = FALSE)
write.csv(lambda_mat, file.path("Exp2_data", "lambda_B_Sample.csv"), row.names = FALSE)

# code is similar with Simulation3, adaptive discretization simulation results
# after import Y_B_Sample
# ---- Load required libraries ----
library(data.table)

# ---- Load observed and true latent data ----
############################ remember to change#################
Y_df <- Y_B_Sample
lambda_true_df <- lambda_B_Sample

TT <- nrow(Y_df)
n_sim <- ncol(Y_df)

# ---- Initialize result containers ----
estimation_results <- list()
Xt_hat_matrix <- matrix(NA, nrow = TT, ncol = n_sim)
lambda_hat_matrix <- matrix(NA, nrow = TT, ncol = n_sim)

# ---- Loop over each simulation column ----
for (i in 1:n_sim) {
  Y <- Y_df[[i]]
  lambda_true <- lambda_true_df[[i]]
  
  # Fit model using local decoding (with error handling)
  results <- tryCatch(
    {
      local_decoding_auto(Y, w)
    },
    error = function(e) {
      warning(sprintf("Error in sim %d of file %s: %s", i, file_name, e$message))
      NULL
    }
  )
  
  if (is.null(results)) {
    init_points <- initial_par(Y, w)
    
    Xt_hat_matrix[, i] <- NA
    lambda_hat_matrix[, i] <- NA
    
    result_row <- data.frame(
      init_phi = tanh(init_points[2]),
      init_sigma = exp(init_points[3]),
      init_beta = exp(init_points[1]),
      est_phi = NA,
      est_sigma = NA,
      est_beta = NA,
      cmin = NA,
      cmax = NA,
      m = NA,
      llk = NA,
      mse = NA,
      nmse = NA
    )
  } else {
    # Extract Xt_hat and lambda_hat
    Xt_hat_matrix[, i] <- results$Xt_hat
    lambda_hat_matrix[, i] <- results$lambda_hat
    
    # Compute mse for lambda_hat
    mse <- mean((results$lambda_hat - lambda_true)^2)
    nmse <- sum((results$lambda_hat - lambda_true)^2) / sum(lambda_true^2)
    
    result_row <- data.frame(
      sim_id = paste0("sim", i),
      init_phi = results$init_points[2],
      init_sigma = results$init_points[3],
      init_beta = results$init_points[1],
      est_phi = results$para_est$phi_hat,
      est_sigma = results$para_est$sigma_hat,
      est_beta = results$para_est$beta_hat,
      cmin = results$cmin,
      cmax = results$cmax,
      m = results$m,
      llk = results$llk,
      mse = mse,
      nmse = nmse
    )
  }
  
  # Store results for each simulation
  estimation_results[[i]] <- result_row
}

# ---- Save outputs ----
est_df <- rbindlist(estimation_results, fill = TRUE)
est_filename <- file.path("Exp2_data", "results.csv")
fwrite(est_df, est_filename)

colnames(Xt_hat_matrix) <- paste0("sim", 1:n_sim)
xhat_filename <- file.path("Exp2_data", "X_hat.csv")
fwrite(as.data.frame(Xt_hat_matrix), xhat_filename)

colnames(lambda_hat_matrix) <- paste0("sim", 1:n_sim)
lhat_filename <- file.path("Exp2_data", "lambda_hat.csv")
fwrite(as.data.frame(lambda_hat_matrix), lhat_filename)

# after importing Exp1_results, if run at once, don't need next line
est_df <- results
beta_est_B <- est_df$est_beta
phi_est_B <- est_df$est_phi
sigma_est_B <- est_df$est_sigma

Theta_mat <- cbind(beta_est_B, phi_est_B, sigma_est_B)
VarCov_hat <- cov(Theta_mat)

# Function to compute equal-tailed 90% CI
boot_ci <- function(x, level = 0.90) {
  alpha <- (1 - level) / 2
  quantile(x, probs = c(alpha, 1 - alpha), na.rm = TRUE)
}

# Compute CIs
beta_ci  <- boot_ci(beta_est_B, level = 0.90)
phi_ci   <- boot_ci(phi_est_B, level = 0.90)
sigma_ci <- boot_ci(sigma_est_B, level = 0.90)

# Print results
beta_ci
phi_ci
sigma_ci
sqrt(diag(VarCov_hat))
mean(results$mse)
table(results$m)
mean(results$llk)


