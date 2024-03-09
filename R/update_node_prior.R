Rcpp::sourceCpp(here::here('src','update_sigma.cpp'))

# Update terminal node prior mean
update_mu_mu <- function (forest, hypers) {
  mu_mu_var_pre <- 100^2
  mu_mu_mean_pre <- 0
  mu <- forest$value[forest$var == -1]
  n_mu <- length(mu)
  mu_mu_var_post <- 1 / ((1 / mu_mu_var_pre) + (n_mu / (hypers$sigma_mu^2)))
  mu_mu_mean_post <- mu_mu_var_post * ((mu_mu_mean_pre / mu_mu_var_pre) + (sum(mu) / (hypers$sigma_mu^2)))
  mu_mu <- stats::rnorm(1, mu_mu_mean_post, sqrt(mu_mu_var_post))

  return(mu_mu)
}

# Update terminal node prior variance
update_sigma_mu <- function (forest, hypers) {

  # # Method 1: Conjugate IG
  # shape <- 0.1 + (sum(forest$var == -1) / 2)
  # scale <- 0.1 + (sum((forest$value[forest$var == -1] - hypers$mu_mu)^2) / 2)
  # sigma_mu <- sqrt(1 / stats::rgamma(1, shape = shape, rate = scale))

  # Method 2: Half-Cauchy
  sigma_mu <- update_sigma(r = forest$value[forest$var == -1] - hypers$mu_mu,
                           sigma_scale = hypers$k / sqrt(hypers$num_tree),
                           sigma_old = hypers$sigma_mu)

  # # Method 3: Marginal Half-Cauchy
  # c <- hypers$k / sqrt(hypers$num_tree)
  # sigma_0_inv <- stats::rgamma(1, shape = 1, rate = (1 / c^2) + (1 / hypers$sigma_mu^2))
  # shape <- (sum(forest$var == -1) + 1) / 2
  # scale <- sum((forest$value[forest$var == -1] - hypers$mu_mu)^2) / 2 + sigma_0_inv
  # sigma_mu <- sqrt(1 / stats::rgamma(1, shape = shape, rate = scale))


  # # Method 5: Inverse Gamma (k)
  # shape <- hypers$k2_shape + (sum(forest$var == -1) / 2)
  # scale <- hypers$k2_scale + (sum((forest$value[forest$var == -1] - hypers$mu_mu)^2) * hypers$num_tree / 2)
  # k2 <- 1 / stats::rgamma(1, shape = shape, rate = scale)
  # sigma_mu <- sqrt(k2 / hypers$num_tree)

  # # Method 6: Marginal half-Cauchy (k)
  # c <- 1
  # k <- hypers$sigma_mu * sqrt(hypers$num_tree)
  # k_0_inv <- stats::rgamma(1, shape = 1, rate = (1 / c^2) + (1 / k^2))
  # shape <- (sum(forest$var == -1) + 1) / 2
  # scale <- (sum((forest$value[forest$var == -1] - hypers$mu_mu)^2) * hypers$num_tree / 2) + k_0_inv
  # k2 <- 1 / stats::rgamma(1, shape = shape, rate = scale)
  # sigma_mu <- sqrt(k2 / hypers$num_tree)

  # # Method 7: Grid search (k)
  # k_grid <- seq(0.01, 3, 0.01)
  # k_post <- sapply(k_grid, \(k) sum(stats::dnorm(forest$value[forest$var == -1], hypers$mu_mu, k / sqrt(hypers$num_tree), log = TRUE)))
  # k <- sample(k_grid, size = 1, prob = exp(k_post - max(k_post)))
  # sigma_mu <- k / sqrt(hypers$num_tree)

  return (sigma_mu)
}



# Update terminal node prior variance using horseshoe prior
update_tau <- function (forest, hypers) {

  scales <- 1 + (1 / hypers$tau^2)
  tau_0_inv <- stats::rgamma(hypers$num_tree, shape = 1, rate = scales)

  tau <- numeric(hypers$num_tree)
  for (t in 1:hypers$num_tree) {
    shape <- (sum(forest$tree == t & forest$var == -1) + 1) / 2
    scale <- sum(((forest$value[forest$tree == t & forest$var == -1] - hypers$mu_mu)^2) / 2 / (hypers$omega^2)) + tau_0_inv[t]
    tau[t] <- sqrt(1 / stats::rgamma(1, shape = shape, rate = scale))
  }

  return (tau)

}
update_omega <- function (forest, hypers) {

  # Draw omega_0 from its full conditional
  scale <- (1 / hypers$sigma_mu_scale^2) + (1 / hypers$omega^2)
  omega_0_inv <- stats::rgamma(1, shape = 1, rate = scale)

  # Draw omega from its full conditional
  mu_tau_sq <- 0
  for (t in 1:hypers$num_tree) {
    mu_sq_t <- (forest$value[forest$tree == t & forest$var == -1] - hypers$mu_mu)^2
    mu_tau_sq <- mu_tau_sq + sum(mu_sq_t / (hypers$tau[t]^2))
  }

  shape <- (sum(forest$var == -1) + 1) / 2
  scale <- (mu_tau_sq / 2) + omega_0_inv

  omega <- sqrt(1 / stats::rgamma(1, shape = shape, rate = scale))

  return (omega)
}



