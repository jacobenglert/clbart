#' Fits a Bayesian conditional logistic regression model
#'
#' Runs a conventional conditional logistic regression model (similar to
#' \code{clogit} in the \code{survival} package) for data which have been
#' processed using the standard (time-stratified) case-crossover approach.
#' @param x A matrix of covariates to include in the linear predictor component.
#' @param y A vector of outcomes (binary).
#' @param strata An integer vector indicating the strata which each observation
#' belongs to. Prior to calling \code{bayes_clr}, all data components must be
#' ordered by strata, and strata should be in increasing order.
#' @param num_burn Number of burn-in iterations for the chain.
#' @param num_thin Post burn-in thinning interval for the chain.
#' @param num_save Number of posterior samples to save; in total, \code{num_burn + num_save * num_thin} iterations are run.
#' @param sigma_beta Proposal distribution variance for fixed effect coefficients.
#' @param update_sigma_beta If \code{TRUE}, \code{sigma_beta} is updated to encourage optimal M-H acceptance rates for the vector of fixed coefficients.
#' @param seed Random seed for reproducibility.
#' @return Returns a list with the following components:
#' - `beta`: \code{num_save} x \code{ncol(x)} matrix of regression coefficients.
#' - `loglik`: vector of length \code{num_save} of the total log-likelihood.
#' - `time`: total time taken for the model to run (in minutes).
#' - `WAIC`: model WAIC (Widely Applicable Information Criteria).
#' @export
#'
#' @examples
#' # TBD
bayes_clr <- function (x, y, strata,
                       num_burn = 2500, num_thin = 1, num_save = 2500,
                       update_sigma_beta = TRUE, sigma_beta = 1,
                       seed = 1996) {

  start <- Sys.time()

  set.seed(seed)

  hypers <- list()
  if (update_sigma_beta) {
    if (num_burn == 0) {
      message("Not updating sigma_beta since no burn-in.")
      update_sigma_beta <- FALSE
    }
    else {
      hypers$beta_update_counter <- 0
      hypers$sigma_beta_update_interval <- num_burn %/% 50
    }
  }


  # Process input data
  x <- as.matrix(x)
  if (is.unsorted(strata)) {
    stop('Sort data by strata before calling function.')
  }
  if (max(strata) != length(unique(strata))) {
    stop('Do not skip over values in strata.')
  }
  p_x <- ncol(x)
  n <- length(unique(strata))


  # Initialize model parameters
  init.fit <- clr_cpp(x, y, strata)
  beta <- init.fit$beta[1:p_x]
  hypers$beta_cov <- init.fit$vcov[1:p_x, 1:p_x]

  # Allocate running storage for WAIC calculation
  ll <- ll2 <- ell <- numeric(n)

  # Prepare posterior storage
  if (num_save > 0) {
    post <- list(

      beta = matrix(NA, ncol = p_x, nrow = num_save, dimnames = list(NULL, colnames(x))),
      sigma_beta = numeric(num_save),
      loglik = numeric(num_save),
      call = match.call()
      )

  } else {
    message("num_save = 0: no posterior samples will be saved.")
  }


  # Run MCMC
  pb <- progress::progress_bar$new(
    format = "[:bar] Burning in :current/:total. Total time elapsed: :elapsedfull",
    total = num_burn, clear = FALSE, width = 100)
  for (k in seq_len(num_burn)) {

    # Update confounders
    beta_old <- beta
    beta <- update_beta(beta, x, y, strata, offset = numeric(nrow(x)), hypers)
    if (update_sigma_beta) {
      if (any(beta != beta_old)) {
        hypers$beta_update_counter <- hypers$beta_update_counter + 1
      }
      if (k %% hypers$sigma_beta_update_interval == 0) {
        hypers$sigma_beta <- update_sigma_beta(hypers)
        hypers$beta_update_counter <- 0
      }
    }

    pb$tick()
  }

  pb <- progress::progress_bar$new(
    format = "[:bar] Sampling from posterior :current/:total. Total time elapsed: :elapsedfull",
    total = num_save, clear = FALSE, width = 100)
  for (k in seq_len(num_save)) {
    for (k2 in seq_len(num_thin)) {
      # Update confounders
      beta <- update_beta(beta, x, y, strata, offset = numeric(nrow(x)), hypers)
    }

    # Store Results
    post$beta[k,]           <- beta
    post$sigma_beta[k]      <- hypers$sigma_beta

    # WAIC calculation
    ll_curr <- clr_loglik_cpp(beta, x, y, strata, offset = numeric(nrow(x)), total = FALSE)
    ll <- ll + ll_curr
    ll2 <- ll2 + ll_curr^2
    ell <- ell + exp(ll_curr)

    post$loglik[k] <- sum(as.numeric(ll_curr))

    pb$tick()
  }

  post$pWAIC <- sum(ll2 / num_save - (ll / num_save)^2)
  post$lppd <- sum(log(ell / num_save))
  post$WAIC <- -2 * post$lppd + 2 * post$pWAIC

  post$time <- as.numeric(difftime(Sys.time(), start, units = 'min'))

  return (post)
}

