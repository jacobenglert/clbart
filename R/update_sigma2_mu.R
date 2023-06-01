update_sigma2_mu <- function(forest, alpha = 0.0001, beta = 0.0001){

  n_mu <- sum(unlist(lapply(forest, '[[', 'num_leaves')))
  mu_sq <- sum(unlist(lapply(forest, \(t) lapply(t$nodes, \(n) n$leaf * (n$mu^2)))))

  return(1 / stats::rgamma(1, shape = n_mu / 2 + alpha, rate = mu_sq / 2 + beta))

}
