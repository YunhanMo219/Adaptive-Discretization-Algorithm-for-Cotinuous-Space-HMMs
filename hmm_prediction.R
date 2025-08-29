# prediction function for poisson hidden  markov model (adapted from Zimmerman et al (2024), Zucchini et al (2016))

# Input: para, estimated parameter (beta, phi, sigma);
#        ess_para, essential domain parameter (cmin, cmax, m);
#        Y, observed data
#        x_list, range of x values; n dim
#        h, 
# Output: h * n matrix


poisson_hmm_forecast <- function(para, ess_para, Y, w, x_list, h){
  sigma <- para$sigma_hat
  phi <- para$phi_hat
  beta <- para$beta_hat
  m <- ess_para$m
  grid_vectors <- discrete_hidden_space(ess_para$cmin, ess_para$cmax, m)
  
  left <- grid_vectors$left
  right <- grid_vectors$right
  mid <- grid_vectors$midpoint
  
  n <- length(x_list)
  TT <- length(Y)
  forecast_matrix <- matrix(0,nrow=h,ncol=n)
  
  phi <- min(phi, 0.999) # avoid phi=1
  phi <- max(phi, -0.999) # avoid phi=-1
  sd_x <- sqrt(sigma^2 / (1 - phi^2))
  pnorm_vals <- pnorm(c(left,right[m]), mean=0, sd=sd_x)
  delta0 <- diff(pnorm_vals)
  delta0 <- delta0 / sum(delta0)
  
  Gam <- matrix(0L, nrow=m, ncol=m)
  for (i in 1:m) {
    for (j in 1:m) {
      Gam[i,j] <- pnorm(right[j], mean=phi*mid[i], sd=sigma) - pnorm(left[j], mean=phi*mid[i], sd=sigma)
    }
  }
  
  lscale <- 0
  ph <- delta0
  for (t in 1:TT) {
    prob <- dpois(Y[t], lambda = w*beta*exp(mid))
    PY <- diag(asFinite(prob))
    v <- as.numeric(ph %*% Gam %*% PY)
    
    # log-sum-exp trick
    log_v <- log(v + 1e-300)
    max_log_v <- max(log_v)
    log_u <- max_log_v + log(sum(exp(log_v - max_log_v)))
    lscale <- lscale + log_u
    
    ph <- exp(log_v - log_u)
    if (lscale==Inf |lscale==-Inf){
      warning("numeric overflow or underflow")
      break
    }
  }
  
  for (i in 1:h) {
    ph <- ph %*% Gam
    ph <- ph / sum(ph)
    for (j in 1:m){
      pois_p <- ph[j] * dpois(x_list, lambda = w*beta*exp(mid)[j])
      forecast_matrix[i,] <- forecast_matrix[i,] + asFinite(pois_p)
    }
  }
  
  return(forecast_matrix)
}