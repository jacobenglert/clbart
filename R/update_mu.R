update_mu <- function(tree, y, z, stratum, sc1 = numeric(length(y)),
                      lambda = numeric(length(y)), sigma2_mu){

  # Log-likelihood and gradient functions needed for ars
  ll <- function(m, y, z, stratum, sc1, lambda, sigma2_mu, max_win, na_locs){
    comp_loglik(y, sc = sc1 + z * (m + lambda), stratum, max_win, na_locs) - ((m^2) / (2*sigma2_mu))
  }
  # grd <- function(m, y, z, stratum, sc1, lambda, sigma2_mu, max_win, na_locs){
  #   comp_grd(y, z, sc1 + z * (m + lambda), stratum, max_win, na_locs) - (m / sigma2_mu)
  # }

  ll_V <- Vectorize(FUN = ll, vectorize.args = 'm')
  # grd_V <- Vectorize(FUN = grd, vectorize.args = 'm')

  leaves <- which(lapply(tree$nodes, \(n) as.numeric(n$leaf)) == 1)
  for(n in leaves){
    id <- tree$nodes[[n]]$data_idx
    # tree$nodes[[n]]$mu <- ars::ars(n = 1, f = ll_V, fprima = grd_V,
    #                                y = y[id], z = z[id], stratum = stratum[id],
    #                                sc1 = sc1[id], lambda = lambda[id],
    #                                sigma2_mu = sigma2_mu,
    #                                max_win = tree$nodes[[n]]$max_win,
    #                                na_locs = tree$nodes[[n]]$na_locs)
    tree$nodes[[n]]$mu <- armspp::arms(n_samples = 1, log_pdf = ll_V,
                                       lower = -1000, upper = 1000,
                                       metropolis = FALSE,
                                       arguments = list(y = y[id], z = z[id],
                                                        stratum = stratum[id],
                                                        sc1 = sc1[id],
                                                        lambda = lambda[id],
                                                        sigma2_mu = sigma2_mu,
                                                        max_win = tree$nodes[[n]]$max_win,
                                                        na_locs = tree$nodes[[n]]$na_locs))
    tree$nodes[[n]]$logLik <- comp_loglik(y = y[id],
                                          sc = sc1[id] + z[id] * (tree$nodes[[n]]$mu + lambda[id]),
                                          stratum = stratum[id],
                                          max_win = tree$nodes[[n]]$max_win,
                                          na_locs = tree$nodes[[n]]$na_locs)
  }

  return(tree)
}
