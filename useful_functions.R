# some useful functions in local decoding function(likelihood_parameter_estimation.R)
# function to left, right, and midpoint vectors for the essential domain grids.
discrete_hidden_space <- function(cmin, cmax, m) {
  # Compute interval length
  delta <- (cmax - cmin) / m
  
  left <- seq(cmin, cmax - delta, length.out = m)
  right <- seq(cmin + delta, cmax, length.out = m)
  midpoint <- (left + right) / 2
  
  df <- data.frame(left = left, right = right, midpoint = midpoint)
  return(df)
}


# function to compute initial parameter values for optimization
initial_par <- function(Y, w, c = 1) {
  # Avoid log(0) by adding a small constant
  Z <- log(Y + c)
  
  # Fit AR(1) to approximate latent process
  ar1_fit <- arima(Z, order = c(1, 0, 0), include.mean = TRUE)
  
  # Extract AR(1) coefficient and innovation standard deviation
  phi_hat <- as.numeric(ar1_fit$coef["ar1"])
  sigma_hat <- sqrt(ar1_fit$sigma2)
  
  # Compute approximate beta from mean(Y) and w
  beta_hat <- mean(Y) / w
  
  beta_hat <- max(beta_hat, 1e-3)
  phi_hat <- min(max(phi_hat, -0.99), 0.99)
  sigma_hat <- max(sigma_hat, 1e-3)
  
  # Return transformed values
  result <- c(
    log_beta = log(beta_hat),
    atanh_phi = atanh(phi_hat),
    log_sigma = log(sigma_hat)
  )
  
  return(result)
}

## functions for log likelihood
asFinite <- function(x) {
  if (any(nifi <- !is.finite(x)))
    x[nifi] <- sign(x[nifi]) * .Machine$double.xmax
  x
}


# function to compute log likelihood of parameters
llk_1d <- function(lbeta, tphi, lsigma, Y, grid_vectors, w){
  beta <- exp(lbeta)
  phi <- tanh(tphi)
  sigma <- exp(lsigma)
  left <- grid_vectors$left
  right <- grid_vectors$right
  mid <- grid_vectors$midpoint
  
  m <- nrow(grid_vectors) # number of grids
  TT <- length(Y)
  
  phi <- min(phi, 0.999) # avoid phi=1
  phi <- max(phi, -0.999) # avoid phi=-1
  sd_x <- sqrt(sigma^2 / (1 - phi^2))
  pnorm_vals <- pnorm(c(left,right[m]), mean=0, sd=sd_x)
  delta0 <- diff(pnorm_vals)
  
  # compute t.p.m
  Gam <- matrix(0L, nrow=m, ncol=m)
  for (i in 1:m) {
    for (j in 1:m) {
      Gam[i,j] <- pnorm(right[j], mean=phi*mid[i], sd=sigma) - pnorm(left[j], mean=phi*mid[i], sd=sigma)
    }
  }
  
  # forward algorithm
  ph <- delta0
  lscale <- 0
  for (t in 1:TT) {
    prob <- dpois(Y[t], lambda = w*beta*exp(mid))
    PY <- diag(asFinite(prob))
    v <- as.numeric(ph %*% Gam %*% PY)
    
    # log-sum-exp trick
    log_v <- log(v + 1e-300)
    max_log_v <- max(log_v)
    log_u <- max_log_v + log(sum(exp(log_v - max_log_v)))
    
    # u <- sum(v)
    # lscale <- lscale + log(asFinite(u))
    lscale <- lscale + log_u
    # ph <- v/u
    
    ph <- exp(log_v - log_u)
    if (lscale==Inf |lscale==-Inf){
      warning("numeric overflow or underflow")
      break
    }
  }
  
  return(asFinite(lscale))
}

# function to compute equal-tailed interval
compute_eti_intervals <- function(stateprobs, states, alpha = 0.95) {
  
  K <- nrow(stateprobs)
  TT <- ncol(stateprobs)
  intervals <- matrix(NA, nrow = TT, ncol = 2)
  
  for (t in 1:TT) {
    probs <- stateprobs[, t]
    probs <- probs / sum(probs)  # 归一化
    cdf <- cumsum(probs)
    
    lower_idx <- which(cdf >= (1 - alpha)/2)[1]
    upper_idx <- which(cdf >= 1 - (1 - alpha)/2)[1]
    
    intervals[t, ] <- c(states[lower_idx], states[upper_idx])
  }
  
  colnames(intervals) <- c("lower", "upper")
  intervals <- as.data.frame(intervals)
  return(intervals)
}


# function to do decoding with estimated parameters and current grids
# add equal-tailed interval
local_decoding <- function(beta_hat, phi_hat, sigma_hat, grid_vectors, Y, w){
  
  Y <- as.matrix(Y)
  m <- nrow(grid_vectors)
  
  left <- grid_vectors$left
  right <- grid_vectors$right
  mid <- grid_vectors$midpoint
  
  delta.hat <- rep(0, m)
  for (j in 1:m) {
    delta.hat[j] <- pnorm(left[j], mean=0, sd=sqrt(sigma_hat^2/(1-phi_hat^2))) 
    - pnorm(right[j], mean=0, sd=sqrt(sigma_hat^2/(1-phi_hat^2)))
  }
  delta.hat <- delta.hat/sum(delta.hat)
  
  Gam.hat <- matrix(0L, nrow=m, ncol=m)
  for (i in 1:m) {
    for (j in 1:m) {
      Gam.hat[i,j] <- pnorm(left[j], mean=phi_hat*mid[i], sd=sigma_hat) 
      - pnorm(right[j], mean=phi_hat*mid[i], sd=sigma_hat)
    }
  }
  Gam.hat <- Gam.hat/rowSums(Gam.hat)
  
  PY_hat <- function(y) {asFinite(dpois(y, lambda = w*beta_hat*exp(mid)))}
  
  
  decoding_results <- local_dec(init=delta.hat, tpm=Gam.hat, data=t(Y), Fx=PY_hat)
  Xt_hat <- mid[decoding_results$preds] # predicted Xt's
  
  stateprob <- decoding_results$probs
  et_interval <- compute_eti_intervals(stateprob, mid, alpha = 0.95)
  
  return(list(Xt_hat=Xt_hat, lower_eti = et_interval$lower, upper_eti = et_interval$upper))
}






