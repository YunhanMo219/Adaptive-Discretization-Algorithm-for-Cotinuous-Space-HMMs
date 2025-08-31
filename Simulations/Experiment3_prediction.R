# Experiment 3 aims at compare the out-of-sample performance of adaptive method 
# and fixed method

# generate data, the setting is similar to earthquake data
# set of the simulation
beta <- 100
phi <- 0.97
sigma <- 0.17
TT <- 125
n_sim <- 100
# b <- 3 * sigma / sqrt(1 - phi^2)
b <- 5
m <- 50

code_dir <- "E:/study in the UK/research project/My Code/Simulation2/Simulation2"
source(file.path(code_dir,"data_simulation.R"))
source(file.path(code_dir,"HMM_Functions.R")) # from Zimmerman et al (2024)
source(file.path(code_dir,"useful_functions.R")) # adapted from Zimmerman et al (2024)
source(file.path(code_dir,"auto_functions.R"))
base_dir2 <- "E:/study in the UK/research project/My Code/Simulation3"
source(file.path(base_dir2, "fixed discretization function.R"))
base_dir3 <- "E:/study in the UK/research project/My Code/Earthquake"
source(file.path(base_dir3, "hmm_prediction.R"))


# Create output folder to store experiment 1 data
dir.create("Exp3_data", showWarnings = FALSE)

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

write.csv(X_mat, file.path("Exp3_data", "X_Sample.csv"), row.names = FALSE)
write.csv(Y_mat, file.path("Exp3_data", "Y_Sample.csv"), row.names = FALSE)
write.csv(lambda_mat, file.path("Exp3_data", "lambda_Sample.csv"), row.names = FALSE)


######### run code from here #########

window_size <- 100
forecast_horizon <- 1
margin <- 40
w <- 1

n_steps <- nrow(Y_Sample) - window_size  # 25
n_series <- ncol(Y_Sample)               # 100

yhat_mat <- matrix(NA, n_steps, n_series)
L95_mat  <- matrix(NA, n_steps, n_series)
U95_mat  <- matrix(NA, n_steps, n_series)
COV_mat  <- matrix(NA, n_steps, n_series)
resid_mat <- matrix(NA, n_steps, n_series)






###########################################################################

for (j in 1:n_series) {
  Y_col <- Y_Sample[[j]]
  for (i in 1:n_steps) {
    Y_window <- Y_col[i:(i + window_size - 1)]
    
    xmin <- max(0, min(Y_window) - margin)
    xmax <- max(Y_window) + margin
    x_list <- xmin:xmax
    
    # res <- local_decoding_fixed(Y = Y_window, w = 1, b = b, m = m) # fixed
    res <- local_decoding_auto(Y = Y_window, w = w) # adaptive
    para <- res$para_est
    ess_para <- list(cmin = res$cmin, cmax = res$cmax, m = res$m)
    
    pmf <- as.numeric(poisson_hmm_forecast(
      para, ess_para, Y_window, w, x_list, forecast_horizon
    )[1, ])
    pmf <- pmf / sum(pmf)
    cat("get prediction!!!!!!! \n")
    
    yhat_mat[i, j] <- x_list[which.max(pmf)]
    
    cdf <- cumsum(pmf)
    low_idx <- which(cdf >= 0.025)[1]
    up_idx  <- which(cdf >= 0.975)[1]
    
    L95_mat[i, j] <- ifelse(is.na(low_idx), xmin, x_list[low_idx])
    U95_mat[i, j] <- ifelse(is.na(up_idx),  xmax, x_list[up_idx])
    
    y_true <- Y_col[i + window_size]
    COV_mat[i, j] <- as.numeric(y_true >= L95_mat[i, j] & y_true <= U95_mat[i, j])
    
    idx_true <- which(x_list == y_true)
    if (length(idx_true) == 0) {
      prob <- ifelse(y_true < xmin, 0, 1)
    } else {
      prob <- cdf[idx_true]
    }
    resid_mat[i, j] <- qnorm(prob)
    cat(j,"th series",i,"th rolling is finished!!!!!!!!!!!!!!!!!! \n")
  }
}

out_dir <- "Exp3_results"
# write.csv(yhat_mat,  file.path(out_dir, "sim_yhat_f_b5_m50.csv"),    row.names = FALSE)
# write.csv(L95_mat,   file.path(out_dir, "sim_Lower95_f_b5_m50.csv"), row.names = FALSE)
# write.csv(U95_mat,   file.path(out_dir, "sim_Upper95_f_b5_m50.csv"), row.names = FALSE)
# write.csv(COV_mat,   file.path(out_dir, "sim_Coverage_f_b5_m50.csv"), row.names = FALSE)
# write.csv(resid_mat, file.path(out_dir, "sim_resid_f_b5_m50.csv"),   row.names = FALSE)

write.csv(yhat_mat,  file.path(out_dir, "sim_yhat_adaptive1.csv"),    row.names = FALSE)
write.csv(L95_mat,   file.path(out_dir, "sim_Lower95_adaptive1.csv"), row.names = FALSE)
write.csv(U95_mat,   file.path(out_dir, "sim_Upper95_adaptive1.csv"), row.names = FALSE)
write.csv(COV_mat,   file.path(out_dir, "sim_Coverage_adaptive1.csv"), row.names = FALSE)
write.csv(resid_mat, file.path(out_dir, "sim_resid_adaptive1.csv"),   row.names = FALSE)


# MSE 
true_mat <- Y_Sample[(window_size+1):nrow(Y_Sample), ]
mse_vec <- colMeans((yhat_mat - true_mat)^2, na.rm=TRUE)
write.csv(data.frame(MSE = mse_vec), file.path(out_dir, "sim_MSE_adaptive1.csv"), row.names = FALSE)











