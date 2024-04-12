#' Create a list of hyperparameter values
#'
#' @param num_tree Number of trees to include in ensemble.
#' @param k Numerator in scale for half-Cauchy prior in terminal node prior standard deviation.
#' @param base Parameter penalizing new nodes in the branching process prior.
#' @param power Parameter penalizing tree depth in the branching process prior.
#' @param mu_mu Mean parameter in terminal node prior.
#' @param sigma_mu Standard deviation parameter in terminal node prior.
#' @param k2_shape Shape parameter for inverse-gamma update of terminal node prior variance (DO NOT USE).
#' @param k2_scale Scale parameter for inverse-gamma update of terminal node prior variance (DO NOT USE).
#' @param s Vector of splitting probabilities for predictors.
#' @param alpha Positive constant controlling the sparsity level.
#' @param alpha_scale Scale of the prior for \code{alpha}; if not provided, defaults to the number of predictors.
#' @param alpha_shape1 Shape parameter for prior on \code{alpha}; if not provided, defaults to 0.5.
#' @param alpha_shape2 Shape parameter for prior on \code{alpha}; if not provided, defaults to 1.0.
#' @param sigma_beta Proposal distribution variance for fixed effect coefficients.
#'
#' @return Returns a list containing the function arguments.
#' @export
Hypers <- function (num_tree = 50, k = 0.5, base = 0.95, power = 2,
                    mu_mu = 0, sigma_mu = NULL,
                    k2_shape = 1.0, k2_scale = 0.25,
                    # omega = NULL, tau = NULL,
                    s = NULL,
                    alpha = 1, alpha_scale = NULL,
                    alpha_shape1 = 0.5, alpha_shape2 = 1,
                    sigma_beta = 1) {

  # Starting value for sigma_mu
  if (is.null(sigma_mu)) sigma_mu <- k / sqrt(num_tree)

  # Ensure probabilities sum to 1
  if (!is.null(s)) s <- s / sum(s)

  out <- list()
  out$num_tree      <- num_tree
  out$k             <- k
  out$base          <- base
  out$power         <- power
  out$mu_mu         <- mu_mu
  out$sigma_mu      <- sigma_mu
  out$s             <- s
  out$alpha         <- alpha
  out$alpha_scale   <- alpha_scale
  out$alpha_shape1  <- alpha_shape1
  out$alpha_shape2  <- alpha_shape2
  out$sigma_beta    <- sigma_beta
  out$k2_shape      <- k2_shape
  out$k2_scale      <- k2_scale

  # # Horseshoe
  # tau_0_inv <- stats::rgamma(num_tree, 1, 2)
  # out$tau <- sqrt(1 / stats::rgamma(num_tree, shape = (num_tree + 1) / 2, rate = tau_0_inv))
  # out$sigma_mu_scale <- k / sqrt(num_tree) # scale of marginal prior
  # out$omega <- k / sqrt(num_tree)
  # out$sigma_mu <- out$omega * out$tau

  return (out)

}

#' MCMC options for clbart
#'
#' @param num_burn Number of burn-in iterations for the chain.
#' @param num_thin Post burn-in thinning interval for the chain.
#' @param num_save Number of posterior samples to save; in total, \code{num_burn + num_save * num_thin} iterations are run.
#' @param num_chains Number of chains to run; values greater than 1 only work with mc_clbart.
#' @param update_sigma_mu If \code{TRUE}, \code{sigma_mu} is updated, with a half-Cauchy prior.
#' @param update_mu_mu If \code{TRUE}, \code{mu_mu} is updated, with a noninformative conjugate normal prior (typically not recommended).
#' @param update_s If \code{TRUE}, \code{s} is updated using the Dirichlet prior \eqn{s \sim D(\alpha / P, \ldots, \alpha / P)} where \eqn{P} is the number of covariates used in the BART component.
#' @param update_alpha If \code{TRUE}, \code{alpha} is updated using a scaled beta prime prior.
#' @param update_sigma_beta If \code{TRUE}, \code{sigma_beta} is updated to encourage optimal M-H acceptance rates for the vector of fixed coefficients.
#' @param store_lambda If \code{TRUE}, values of the BART predictor for each individual are stored. In practice, this can be heavy on system memory.
#' @param seed Random seed for reproducibility.
#'
#' @return Returns a list containing the function arguments
#' @export
Opts <- function (num_burn = 2500, num_thin = 1, num_save = 2500,
                  num_chains = 1,
                  update_sigma_mu = TRUE,
                  update_mu_mu = FALSE,
                  update_s = FALSE,
                  update_alpha = FALSE,
                  update_sigma_beta = TRUE,
                  store_lambda = FALSE,
                  seed = 1996) {

  if (num_save < 0 | !(num_save %% 1 == 0)) stop('num_save must be a non-negative integer.')
  if (num_burn < 0 | !(num_burn %% 1 == 0)) stop('num_burn must be a non-negative integer.')
  if (num_thin < 1 | !(num_thin %% 1 == 0)) stop('num_thin must be a positive integer.')

  out <- list()
  out$num_burn <- num_burn
  out$num_thin <- num_thin
  out$num_save <- num_save
  out$num_chains <- num_chains
  out$update_mu_mu <- update_mu_mu
  out$update_sigma_mu <- update_sigma_mu
  out$update_s <- update_s
  out$update_alpha <- update_alpha
  out$update_sigma_beta <- update_sigma_beta
  out$store_lambda <- store_lambda
  out$seed <- seed


  return (out)
}


#' Fits the clbart model
#'
#' Runs the CL-BART model for data which have been processed using the standard
#' (time-stratified) case-crossover approach.
#'
#' @param w A matrix of covariates to include in the BART component.
#' @param x A matrix of covariates to include in the linear predictor component.
#' @param y A vector of outcomes (binary).
#' @param z A vector of primary exposures. This may be continuous or binary, but
#' must vary within each strata or else errors will occur.
#' @param strata An integer vector indicating the strata which each observation
#' belongs to. Prior to calling \code{clbart}, all data components must be
#' ordered by strata, and strata should be in increasing order.
#' @param hypers A list of hyperparameter values obtained from the \code{Hypers} function.
#' @param opts A list of MCMC options obtained from the \code{Opts} function.
#'
#' @return Returns a list with the following components:
#' @export
#'
#' @examples
#' # TBD
clbart <- function (w, x, y, z, strata, hypers, opts) {

  set.seed(opts$seed)

  if (is.null(hypers$s)) hypers$s <- rep(1 / ncol(w), ncol(w))
  if (opts$update_s) hypers$logs <- log(hypers$s)
  if (opts$update_alpha & is.null(hypers$alpha_scale)) hypers$alpha_scale <- ncol(w)

  if (opts$update_sigma_beta) {
    if (opts$num_burn == 0) {
      message("Not updating sigma_beta since no burn-in.")
      opts$update_sigma_beta <- FALSE
    }
    else {
      hypers$beta_update_counter <- 0
      hypers$sigma_beta_update_interval <- opts$num_burn %/% 50
    }
  }


  # Process input data
  w <- as.matrix(w)
  x <- as.matrix(x)
  z <- matrix(z, ncol = 1)
  if (is.unsorted(strata)) {
    stop('Sort data by strata before calling function.')
  }
  if (max(strata) != length(unique(strata))) {
    stop('Do not skip over values in strata.')
  }
  p_x <- ncol(x)
  p_w <- ncol(w)
  n <- length(unique(strata))


  # Initialize model parameters
  # Confounders
  init.fit <- clr_cpp(cbind(x, z), y, strata)
  beta <- init.fit$beta[1:p_x]
  hypers$beta_cov <- init.fit$vcov[1:p_x, 1:p_x]
  # hypers$mu_mu <- init.fit$beta[p_x + 1]

  # Forest
  forest <- make_forest(w, y, z, strata, offset = x %*% beta, hypers)
  lambda <- rep(sum(forest$value), nrow(w))

  # Allocate running storage for WAIC calculation
  ll <- ll2 <- ell <- numeric(n)

  # Prepare posterior storage
  if (opts$num_save > 0) {
    post <- list(

      beta = matrix(NA, ncol = p_x, nrow = opts$num_save, dimnames = list(NULL, colnames(x))),
      sigma_beta = numeric(opts$num_save),
      sigma_mu = numeric(opts$num_save),
      mu_mu = numeric(opts$num_save),
      alpha = numeric(opts$num_save),

      avg_num_nodes = numeric(opts$num_save),
      avg_tree_depth = numeric(opts$num_save),
      split_counts = matrix(NA, ncol = p_w, nrow = opts$num_save, dimnames = list(NULL, colnames(w))),
      split_probs = matrix(NA, ncol = p_w, nrow = opts$num_save, dimnames = list(NULL, colnames(w))),

      loglik = numeric(opts$num_save),

      # omega = numeric(opts$num_save),
      # tau = matrix(NA, ncol = hypers$num_tree, nrow = opts$num_save),

      lambda_mean_overall = numeric(opts$num_save),
      lambda_est = numeric(nrow(w)),

      forest = list(),

      call = match.call()
      )

    if (opts$store_lambda) post$lambda <- matrix(NA, ncol = nrow(w), nrow = opts$num_save)

  } else {
    message("num_save = 0: no posterior samples will be saved.")
  }


  # Run MCMC
  pb <- progress::progress_bar$new(
    format = "[:bar] Burning in :current/:total. Total time elapsed: :elapsedfull",
    total = opts$num_burn, clear = FALSE, width = 100)
  for (k in seq_len(opts$num_burn)) {

    # Update confounders
    beta_old <- beta
    beta <- update_beta(beta, x, y, strata, z * lambda, hypers)
    if (opts$update_sigma_beta) {
      if (any(beta != beta_old)) {
        hypers$beta_update_counter <- hypers$beta_update_counter + 1
      }
      if (k %% hypers$sigma_beta_update_interval == 0) {
        hypers$sigma_beta <- update_sigma_beta(hypers)
        hypers$beta_update_counter <- 0
      }
    }
    x_beta <- x %*% beta

    # Update forest
    forest <- update_forest(forest, w, y, z, strata, x_beta, hypers)
    lambda <- get_forest_predictions(forest, w)

    # Update other parameters
    if (k > opts$num_burn / 2) {
      if (opts$update_s) {
        hypers$logs <- update_s(forest, hypers)
        hypers$s <- logs_to_s(hypers$logs)
      }
      if (opts$update_alpha) hypers$alpha <- update_alpha(hypers)
    }

    if (opts$update_mu_mu) hypers$mu_mu <- update_mu_mu(forest, hypers)
    if (opts$update_sigma_mu) hypers$sigma_mu <- update_sigma_mu(forest, hypers)
    # hypers$tau <- update_tau(forest, hypers)
    # hypers$omega <- update_omega(forest, hypers)

    pb$tick()
  }

  pb <- progress::progress_bar$new(
    format = "[:bar] Sampling from posterior :current/:total. Total time elapsed: :elapsedfull",
    total = opts$num_save, clear = FALSE, width = 100)
  for (k in seq_len(opts$num_save)) {
    for (k2 in seq_len(opts$num_thin)) {
      # Update confounders
      beta <- update_beta(beta, x, y, strata, z * lambda, hypers)
      x_beta <- x %*% beta

      # Update forest
      forest <- update_forest(forest, w, y, z, strata, x_beta, hypers)
      lambda <- get_forest_predictions(forest, w)

      # Update other parameters
      if (opts$update_s) {
        hypers$logs <- update_s(forest, hypers)
        hypers$s <- logs_to_s(hypers$logs)
      }
      if (opts$update_alpha) hypers$alpha <- update_alpha(hypers)
      if (opts$update_mu_mu) hypers$mu_mu <- update_mu_mu(forest, hypers)
      if (opts$update_sigma_mu) hypers$sigma_mu <- update_sigma_mu(forest, hypers)
      # hypers$tau <- update_tau(forest, hypers)
      # hypers$omega <- update_omega(forest, hypers)

    }

    # Store Results
    post$alpha[k]           <- hypers$alpha
    post$beta[k,]           <- beta

    post$avg_num_nodes[k]   <- mean(table(forest$tree))
    post$avg_tree_depth[k]  <- mean(tapply(forest$node, forest$tree, \(x) get_node_depth(max(x))))
    post$split_counts[k,]   <- sapply(1:p_w, \(j) sum(forest$var == j))
    post$split_probs[k,]    <- hypers$s

    post$mu_mu[k]           <- hypers$mu_mu
    post$sigma_mu[k]        <- hypers$sigma_mu
    post$sigma_beta[k]      <- hypers$sigma_beta

    # post$omega[k] <- hypers$omega
    # post$tau[k,] <- hypers$tau

    post$forest[[k]] <- cbind(sample = rep.int(k, nrow(forest)),
                              forest[,c('tree','node','n','var','value','loglik')])

    if (opts$store_lambda) post$lambda[k,] <- lambda
    post$lambda_mean_overall[k] <- mean(lambda)
    post$lambda_est <- post$lambda_est + (lambda / opts$num_save)

    # WAIC calculation
    ll_curr <- clr_loglik_cpp(0, z, y, strata, offset = x_beta + z * lambda, total = FALSE)
    ll <- ll + ll_curr
    ll2 <- ll2 + ll_curr^2
    ell <- ell + exp(ll_curr)

    post$loglik[k] <- sum(as.numeric(ll_curr))

    pb$tick()
  }

  post$forest <- do.call(rbind, post$forest)
  post$pWAIC <- sum(ll2 / opts$num_save - (ll / opts$num_save)^2)
  post$lppd <- sum(log(ell / opts$num_save))
  post$WAIC <- -2 * post$lppd + 2 * post$pWAIC

  return (post)
}

