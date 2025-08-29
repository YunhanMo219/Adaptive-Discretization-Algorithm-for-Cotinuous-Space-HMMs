# set of the simulation
beta <- 20
phi <- 0.9
sigma <- 0.2
# sigma <- 0.4
TT <- 1000

#b_H<- 3 * sigma / sqrt(1 - phi^2)
b <- 1
m <- 100
n_sim <- 200

# source functinos
simu_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"
source(file.path(simu_dir, "data_simulation.R"))

code_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"
source(file.path(code_dir,"data_simulation.R"))
source(file.path(code_dir,"HMM_Functions.R")) # from Zimmerman et al (2024)
source(file.path(code_dir,"useful_functions.R")) # adapted from Zimmerman et al (2024)
source(file.path(code_dir,"auto_functions.R"))

code_dir2 <- "E:/study in the UK/research project/My Code/Earthquake"
source(file.path(code_dir2, "hmm_prediction.R"))

# after import Y_B_Sample
# ---- Load required libraries ----
library(data.table)

# ---- Load observed and true latent data ----
############## remember to change data################3
Y_df <- Y_B_Sample
lambda_true_df <- lambda_B_Sample

TT <- nrow(Y_df)
n_sim <- ncol(Y_df)

# ---- Initialize result containers ----
estimation_results_f1 <- list()
Xt_hat_matrix_f1 <- matrix(NA, nrow = TT, ncol = n_sim)
lambda_hat_matrix_f1 <- matrix(NA, nrow = TT, ncol = n_sim)

# ---- Loop over each simulation column ----
for (i in 1:n_sim) {
  Y <- Y_df[[i]]
  lambda_true <- lambda_true_df[[i]]
  
  # Fit model using local decoding (with error handling)
  results <- tryCatch(
    {
      local_decoding_fixed(Y, w=1, b=b, m=m)
    },
    error = function(e) {
      warning(sprintf("Error in sim %d of file %s: %s", i, file_name, e$message))
      NULL
    }
  )
  
  if (is.null(results)) {
    init_points <- initial_par(Y, w)
    
    Xt_hat_matrix_f1[, i] <- NA
    lambda_hat_matrix_f1[, i] <- NA
    
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
    Xt_hat_matrix_f1[, i] <- results$Xt_hat
    lambda_hat_matrix_f1[, i] <- results$lambda_hat
    
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
  estimation_results_f1[[i]] <- result_row
}

# ---- Save outputs ----
est_df_f1 <- rbindlist(estimation_results_f1, fill = TRUE)
est_filename_f1 <- file.path("Exp1_data", "results_fixed_b1_m100.csv")
fwrite(est_df_f1, est_filename_f1)

colnames(Xt_hat_matrix_f1) <- paste0("sim", 1:n_sim)
xhat_filename_f1 <- file.path("Exp1_data", "X_hat_fixed_b1_m100.csv")
fwrite(as.data.frame(Xt_hat_matrix_f1), xhat_filename_f1)

colnames(lambda_hat_matrix_f1) <- paste0("sim", 1:n_sim)
lhat_filename_f1 <- file.path("Exp1_data", "lambda_hat_fixed_b1_m100.csv")
fwrite(as.data.frame(lambda_hat_matrix_f1), lhat_filename_f1)




# after importing Exp1_results, if run at once, don't need next line
est_df_f1 <- results_fixed_b1_m100
beta_est_B_f1 <- est_df_f1$est_beta
phi_est_B_f1 <- est_df_f1$est_phi
sigma_est_B_f1 <- est_df_f1$est_sigma

Theta_mat_f1 <- cbind(beta_est_B_f1, phi_est_B_f1, sigma_est_B_f1)
VarCov_hat_f1 <- cov(Theta_mat_f1)

# Function to compute equal-tailed 90% CI
boot_ci <- function(x, level = 0.90) {
  alpha <- (1 - level) / 2
  quantile(x, probs = c(alpha, 1 - alpha), na.rm = TRUE)
}

# Compute CIs
beta_ci_f1  <- boot_ci(beta_est_B_f1, level = 0.90)
phi_ci_f1   <- boot_ci(phi_est_B_f1, level = 0.90)
sigma_ci_f1 <- boot_ci(sigma_est_B_f1, level = 0.90)

# Print results
beta_ci_f1
phi_ci_f1
sigma_ci_f1
sqrt(diag(VarCov_hat_f1))
mean(est_df_f1$mse)
mean(est_df_f1$llk)
