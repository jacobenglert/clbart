update_mu <- function (tree, w, y, z, s, offset, hypers) {

  # Identify leaves
  leaves <- tree[tree$var == -1,]

  # Iterate through each leaf
  for (i in 1:nrow(leaves)) {

    # Subset leaf
    l <- leaves[i,]

    # Retrieve observations mapped to l
    l.ind <- get_mapped_observations(tree, l$node, w)
    s1 <- match(s[l.ind], unique(s[l.ind]))
    x1 <- z[l.ind,,drop = FALSE]
    y1 <- y[l.ind]
    o1 <- offset[l.ind]

    # Create vectorized log-density function for l
    log_posterior <- function(value) {
      clr_loglik_cpp(value, x1, y1, s1, o1) - (((value - hypers$mu_mu)^2) / (2 * hypers$sigma_mu^2))
    }
    log_posterior_v <- Vectorize(FUN = log_posterior, vectorize.args = 'value')

    # gradient <- function(value) {
    #   clr_grd_cpp(value, x1, y1, s1, o1) - ((value - hypers$mu_mu) / (hypers$sigma_mu^2))
    # }
    # gradient_v <- Vectorize(FUN = gradient, vectorize.args = 'value')
    #
    # # Sample using ARS
    # l$value <- ars::ars(n = 1, f = log_posterior_v, fprima = gradient_v)

    # Sample using ARS
    l$value <- armspp::arms(n_samples = 1,
                            log_pdf = log_posterior_v,
                            lower = -10, upper = 10,
                            # lower = l$m - 3*l$v, upper = l$m + 3*l$v,
                            metropolis = FALSE)

    # l$value <- slice_sampler(log_posterior_v, start = l$m, steps = 1, w = 6 * l$v, m = 20)

    # Recalculate node log-likelihood
    l$loglik <- as.numeric(clr_loglik_cpp(l$value, x1, y1, s1, o1))

    # Update tree
    tree[tree$node == l$node,] <- l
  }

  return(tree)
}


# slice_sampler <- function(f, start, steps=100, w=1, m=10) {
#   # f: Scalar-valued target distribution function
#   # start: Initial value to start the sampling
#   # steps: Number of MCMC steps to perform
#   # w: Initial width of the slice sampling interval
#   # m: Maximum number of steps outward the interval can be expanded
#
#   samples <- numeric(steps)  # Vector to store samples
#   current_x <- start          # Current sample position
#
#   for (i in 1:steps) {
#     y <- f(current_x) - stats::rexp(1)  # Draw a horizontal slice level (log scale)
#
#     # Define an interval [L, R] around current_x
#     L <- current_x - stats::runif(1) * w
#     R <- L + w
#
#     # Step out the interval until the end points are outside the slice
#     for (j in 1:m) {
#       if (f(L) < y) break
#       L <- L - w
#     }
#     for (j in 1:m) {
#       if (f(R) < y) break
#       R <- R + w
#     }
#
#     # Sample from the interval [L, R] until a valid sample is found
#     repeat {
#       x_proposal <- runif(1, L, R)  # Uniformly propose a new point in the interval
#       if (f(x_proposal) >= y) {     # Check if the proposed point is above the slice
#         current_x <- x_proposal     # Accept the proposal
#         break
#       } else {                      # Otherwise, shrink the interval
#         if (x_proposal < current_x) {
#           L <- x_proposal
#         } else {
#           R <- x_proposal
#         }
#       }
#     }
#
#     samples[i] <- current_x  # Store the sample
#   }
#
#   samples
# }



# mu_range <- seq(l$m - 2*l$v, l$m + 2*l$v, by = .01)
# mu_range <- seq(-.2,0.2,0.01)
# ll_range <- log_posterior_v(mu_range)
# plot(mu_range, (ll_range - min(ll_range)) / (max(ll_range) - min(ll_range)),
#      type = 'l', xlab = 'mu', lwd = 2, col = 'blue')
# mu_samples <- slice_sampler(log_posterior_v, start = l$m, step = 1000, w = 4*l$v)
# hist(mu_samples, breaks = 30, probability = TRUE)
#
#
# mu_samps_ars <- ars::ars(n = 1000, f = log_posterior_v, fprima = gradient_v)
#                     #   x = c(l$m - 2*l$v, l$m, l$m + 2*l$v))
#
# mu_samps_armspp <- armspp::arms(n_samples = 1000,
#                              log_pdf = log_posterior_v,
#                              lower = l$m - 3*l$v, upper = l$m + 3*l$v,
#                              metropolis = FALSE)
#
# mu_samps_armsppm <- armspp::arms(n_samples = 1000,
#                                 log_pdf = log_posterior_v,
#                                 lower = hypers$mu_mu - 3 * hypers$sigma_mu,
#                                 upper = hypers$mu_mu + 3 * hypers$sigma_mu,
#                                 metropolis = TRUE)
#
#
# mu_samps_slice <- slice_sampler(log_posterior_v, start = l$m, step = 1000, w = 6*l$v, m = 10)
#
# par(mfrow = c(2,2))
# hist(mu_samps_ars)
# hist(mu_samps_armspp)
# hist(mu_samps_armsppm)
# hist(mu_samps_slice)
#
# get_node_map(z, y, s, offset, l.ind, hypers)
# mean(mu_samps_armspp)
# mean(mu_samps_armsppm)
# mean(mu_samps_ars)
# mean(mu_samps_slice)
#
#
#
# rbenchmark::benchmark(replications = 100,
#                       'armspp' = armspp::arms(n_samples = 1,
#                                               log_pdf = log_posterior_v,
#                                               lower = l$m - 3*l$v, upper = l$m + 3*l$v,
#                                               metropolis = FALSE),
#                       'armsppm' = armspp::arms(n_samples = 1,
#                                                log_pdf = log_posterior_v,
#                                                lower = l$m - 3*l$v, upper = l$m + 3*l$v,
#                                                metropolis = TRUE),
#                       'slice' = slice_sampler(log_posterior_v, start = l$m, steps = 1, w = 6*l$v, m = 5),
#                       'slice_cpp' = slice_sampler_cpp(log_posterior_v, start = l$m, steps = 1, w = 6*l$v, m = 5))

