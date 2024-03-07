update_beta <- function(beta, x, y, strata, offset, hypers) {
  
  prior_cov <- 1e10 * diag(length(beta))
  
  # Betas (proposal and current)
  beta_curr <- as.numeric(beta)
  beta_prop <- as.numeric(mvtnorm::rmvnorm(1, beta_curr, hypers$sigma_beta^2 * hypers$beta_cov))
  
  # Log-likelihoods (proposal and current)
  ll_curr <- clr_loglik_cpp(beta_curr, x, y, strata, offset)
  ll_prop <- clr_loglik_cpp(beta_prop, x, y, strata, offset)
  
  # Log-priors (proposal and current)
  lp_curr <- mvtnorm::dmvnorm(beta_curr, sigma = prior_cov, log = TRUE)
  lp_prop <- mvtnorm::dmvnorm(beta_prop, sigma = prior_cov, log = TRUE)
  
  # Acceptance ratio
  r <- min(ll_prop + lp_prop - ll_curr - lp_curr, 0)
  
  # Update beta vector (or not)
  if (r > log(stats::runif(1))) {
    beta <- beta_prop
  } else {
    beta <- beta_curr
  }
  
  return(beta)
}

update_sigma_beta <- function (hypers) {

    rate <- hypers$beta_update_counter / hypers$sigma_beta_update_interval
    
    if (rate < 0.001) sigma_beta <- hypers$sigma_beta * sqrt(0.1)
    else if (rate < 0.05) sigma_beta <- hypers$sigma_beta * sqrt(0.5)
    else if (rate < 0.20) sigma_beta <- hypers$sigma_beta * sqrt(0.9)
    else if (rate > 0.50) sigma_beta <- hypers$sigma_beta * sqrt(1.1)
    else if (rate > 0.75) sigma_beta <- hypers$sigma_beta * sqrt(2.0)
    else if (rate > 0.95) sigma_beta <- hypers$sigma_beta * sqrt(10.0)
    else sigma_beta <- hypers$sigma_beta
    
    return(sigma_beta)
}
