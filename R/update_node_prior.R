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

  # Method 1: Half-Cauchy
  sigma_mu <- update_sigma(r = forest$value[forest$var == -1] - hypers$mu_mu,
                           sigma_scale = hypers$k / sqrt(hypers$num_tree),
                           sigma_old = hypers$sigma_mu)

  # # Method 2: Marginal Half-Cauchy
  # c <- hypers$k / sqrt(hypers$num_tree)
  # sigma_0_inv <- stats::rgamma(1, shape = 1, rate = (1 / c^2) + (1 / hypers$sigma_mu^2))
  # shape <- (sum(forest$var == -1) + 1) / 2
  # scale <- sum((forest$value[forest$var == -1] - hypers$mu_mu)^2) / 2 + sigma_0_inv
  # sigma_mu <- sqrt(1 / stats::rgamma(1, shape = shape, rate = scale))

  # # Method 3: Marginal half-Cauchy (k)
  # c <- 1
  # k <- hypers$sigma_mu * sqrt(hypers$num_tree)
  # k_0_inv <- stats::rgamma(1, shape = 1, rate = (1 / c^2) + (1 / k^2))
  # shape <- (sum(forest$var == -1) + 1) / 2
  # scale <- (sum((forest$value[forest$var == -1] - hypers$mu_mu)^2) * hypers$num_tree / 2) + k_0_inv
  # k2 <- 1 / stats::rgamma(1, shape = shape, rate = scale)
  # sigma_mu <- sqrt(k2 / hypers$num_tree)

  return (sigma_mu)
}
