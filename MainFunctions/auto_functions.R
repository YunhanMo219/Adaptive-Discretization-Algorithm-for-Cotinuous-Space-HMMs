# two functions are defined in this file, one to estimate parameters under likelihood framework,
# one to find the enssential domain and m automatically.

# function to estimate parameters under likelihood framework
likelihood_est <- function(Y, grid_vectors, w, par_init = c(log(100), atanh(0.9),log(2))){
  
  # compute the log likelihood based on the current grids
  # par_init <- c(log(0.6), atanh(0.9),log(0.3)) # Transform the parameters to convert the problem into an unconstrained optimization
  
  # estimate parameter by maxing log likelihood using optimization methods
  fit <- optim(
    par = par_init,
    fn = function(p) -llk_1d(lbeta=p[1],
                             tphi=p[2],
                             lsigma=p[3],
                             Y,
                             grid_vectors,
                             w), # llk_1d was defined in useful_functions.R
    method = "BFGS",
    control = list(trace = 1),
    hessian=FALSE
  )
  
  # get estimated parameters
  beta_hat <- exp(fit$par[1])
  phi_hat <- tanh(fit$par[2])
  sigma_hat <- exp(fit$par[3])
  
  # compute llk
  llk_1d <- llk_1d(fit$par[1], fit$par[2], fit$par[3], Y, grid_vectors, w)
  
  return(list(beta_hat = beta_hat, phi_hat = phi_hat, sigma_hat = sigma_hat, llk = llk_1d))
}


#   Likelihood-based decoding with automatic essential domain and m selection.
#   For essential domain, The essential domain is chosen such that the maximum 
# and minimum values of the decoded result do not fall on the boundary grid points; 
# if they do, the domain is expanded accordingly.
#   For m, Increase m until the maximized likelihood stays constant and 
# parameter estimates are stable.

# Input: observed data Y
# Output: essential domain, m, decoded results.


local_decoding_auto <- function(Y, w, max_adjust = 10, max_m_adjust = 10, verbose = TRUE){
  
  # set initial essential domain and m
  m <- 50 # the number of grids in the essential domain
  init_points <- initial_par(Y, w) # get some initial estimation of parameters
  # init_points <- c(log(30), atanh(0.9),log(0.5))
  # b <- 10 * exp(init_points[3]) # the initial essential domain is set as [-initial_sigma, initial_sigma]
  b <- 5
  cmin <- -b
  cmax <- b
  
  # tune essential domain
  adjust_counter <- 0 # for essential domain
  converged <- FALSE
  
  while (!converged && adjust_counter < max_adjust){
    adjust_counter <- adjust_counter + 1
    delta_d <- (cmax - cmin) / m
    expand <- FALSE
    
    # compute the left, right, and center points of each grid
    grid_vectors <- discrete_hidden_space(cmin, cmax, m) # defined in useful_functions.R
    K <- nrow(grid_vectors) # number of grid points
    
    # estimate parameters (using the computed initial points)
    para_est <- likelihood_est(Y, grid_vectors, w, par_init = init_points)
    
    # decode hidden state
    Xt_hat_all <- local_decoding(
      para_est$beta_hat, 
      para_est$phi_hat, 
      para_est$sigma_hat, 
      grid_vectors, 
      Y,
      w)
    
    Xt_hat <- Xt_hat_all$Xt_hat
    
    if (max(Xt_hat) == grid_vectors$midpoint[K]) {
      cmax <- cmax + 0.5
      expand <- TRUE
      if (verbose) cat("Expand cmax to:", cmax, "\n")
    } else {
      cmax_new <- max(Xt_hat) + delta_d
      if (cmax_new != cmax) {
        cmax <- cmax_new
        if (verbose) cat("Shrink cmax to:", cmax, "\n")
      }
    }
    
    if (min(Xt_hat) == grid_vectors$midpoint[1]) {
      cmin <- cmin - 0.5
      expand <- TRUE
      if (verbose) cat("Expand cmin to:", cmin, "\n")
    } else {
      cmin_new <- min(Xt_hat) - delta_d
      if (cmin_new != cmin) {
        cmin <- cmin_new
        if (verbose) cat("Shrink cmin to:", cmin, "\n")
      }
    }
    
    if (!expand) {
      converged <- TRUE
      if (verbose) cat("Essential domain converged.\n")
    }
  }
  
  if (adjust_counter >= max_adjust && !converged) {
    warning("Essential domain did not converge within the maximum number of adjustments.")
  }
  
  # tune m
  m_adjust_counter <- 0
  likelihood_diff <- Inf
  likelihood_threshold <- 1e-4
  
  loglik_prev <- para_est$llk
  
  while (m_adjust_counter < max_m_adjust && likelihood_diff > likelihood_threshold) {
    m_adjust_counter <- m_adjust_counter + 1
    
    m <- m + 50
    if (verbose) cat("Increasing m to:", m, "\n")
    
    grid_vectors <- discrete_hidden_space(cmin, cmax, m)
    
    para_est_new <- likelihood_est(Y, grid_vectors, w, par_init = init_points)
    loglik_new <- para_est_new$llk
    
    likelihood_diff <- abs(loglik_new - loglik_prev) / (abs(loglik_prev) + 1e-10)
    if (verbose) cat("Log-likelihood diff:", likelihood_diff, "\n")
    
    prev_para_est <- para_est
    para_est <- para_est_new
    loglik_prev <- loglik_new
  }
  
  if (m_adjust_counter >= max_m_adjust && likelihood_diff > likelihood_threshold) {
    warning("m adjustment reached the maximum number without stable likelihood.")
  }
  
  if (m_adjust_counter < max_m_adjust) {
    m <- m - 50
    cat("Final m is:", m, "\n")
  }
  
  grid_vectors <- discrete_hidden_space(cmin, cmax, m)
  Xt_hat_all <- local_decoding(
    prev_para_est$beta_hat, 
    prev_para_est$phi_hat, 
    prev_para_est$sigma_hat, 
    grid_vectors, 
    Y,
    w)
  
  Xt_hat <- Xt_hat_all$Xt_hat
  
  init_points[1] <- exp(init_points[1])
  init_points[2] <- tanh(init_points[2])
  init_points[3] <- exp(init_points[3])
  
  # more results
  # make sure that Xt_hat has zero mean
  Xt_hat_mean <- mean(Xt_hat)
  Xt_hat <- Xt_hat - Xt_hat_mean
  para_est$beta_hat <- para_est$beta_hat * exp(Xt_hat_mean)
  cmin <- cmin - Xt_hat_mean
  cmax <- cmax - Xt_hat_mean
  
  Xt_lower_eti <- Xt_hat_all$lower_eti - Xt_hat_mean
  Xt_upper_eti <- Xt_hat_all$upper_eti - Xt_hat_mean
  
  lambda_hat <- w * para_est$beta_hat * exp(Xt_hat)
  lambda_lower_eti <- w * para_est$beta_hat * exp(Xt_lower_eti)
  lambda_upper_eti <- w * para_est$beta_hat * exp(Xt_upper_eti)
  
  results <- list(cmin = cmin, cmax = cmax, m = m, Xt_hat = Xt_hat, 
                  adjust_ed_count = adjust_counter, adjust_m_count = m_adjust_counter, 
                  para_est = para_est, init_points = init_points, llk=loglik_prev,
                  lambda_hat = lambda_hat, lambda_lower_eti = lambda_lower_eti, 
                  lambda_upper_eti = lambda_upper_eti)
  return(results)
}

