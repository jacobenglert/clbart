update_alpha <- function (hypers) {

  p_w <- length(hypers$s)
  # l_s_w <- log(s_w)

  theta_grid <- seq.int(1, 999) / 1000

  # Compute posterior density over theta grid
  ll_theta <- sapply(theta_grid, \(theta) {

    # Reparamaterize in terms of alpha
    alpha <- hypers$alpha_scale * theta / (1 - theta)

    # Dirichlet log-likelihood
    log_lik <- lgamma(alpha) - p_w * lgamma(alpha / p_w) + (alpha / p_w) * sum(hypers$logs)

    # log beta prior density
    log_prior <- stats::dbeta(theta, shape1 = hypers$alpha_shape1, shape2 = hypers$alpha_shape2, log = TRUE)

    # log posterior density
    log_post <- log_lik + log_prior
    return (log_post)
  })

  # Sample theta and convert to alpha
  theta <- sample(theta_grid, 1, prob = exp(ll_theta - log_sum_exp(ll_theta)))
  alpha <- theta * hypers$alpha_scale / (1 - theta)

  return(alpha)
}

update_s <- function (forest, hypers) {

  p_w <- length(hypers$s)

  # Tally observed splits
  u <- sapply(1:p_w, \(j) sum(forest$var == j))

  # # Attempted splits before finding a valid one
  # b_valid_vars <- forest$goodvars[forest$var != -1]
  # Z <- sapply(b_valid_vars, \(b) {
  #
  #   p_w <- length(b)
  #
  #   # If all variables have valid splits, return 0 vector
  #   if (all(unlist(b))) return(numeric(p_w))
  #
  #   # Get splitting probabilities for ineligible predictors only
  #   s_w_bad <- exp(hypers$logs) * !unlist(b)
  #
  #   # Sample a number of attempted splits made
  #   if (round(sum(s_w_bad), 5) == 1) {
  #     # Since the node was split, the splitting probabilties have shifted such
  #     # that an infinite number of latent split attempts would be required now,
  #     # so just add 1 for each invalid split
  #     return(as.numeric(!unlist(b)))
  #   } else {
  #
  #     # Sample the number of latent split attempts prior to selecting a good one
  #     W <- rgeom(1, 1 - sum(s_w_bad))
  #
  #     # If zero latent splits were attempted, return 0 vector
  #     if (W == 0) return(numeric(p_w))
  #
  #     # Sample a distribution of latent split attempts across invalid variables
  #     Z <- rowSums(rmultinom(W, 1, s_w_bad)) #s_w_b))
  #     return(Z)
  #   }
  # })
  # u <- rowSums(Z) + u

  # Sample s_w on the log scale
  l_s_w <- sapply(hypers$alpha / p_w + u, rlgam)

  # Normalize log s_w values
  l_s_w <- l_s_w - log_sum_exp(l_s_w)

  # Convert to probability scale
  # s_w <- exp(l_s_w)

  return(l_s_w)
}
