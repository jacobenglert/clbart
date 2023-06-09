#' Estimate Heterogeneous Effects with Conditional Logistic Regression
#'
#' @description
#' Estimate heterogeneous effects of a single exposure/treatment covariate with
#' conditional logistic regression.
#'
#'
#' @param w A data frame of potential effect moderators.
#' @param x A data frame of effect confounders.
#' @param y A numeric vector of 1's and 0's corresponding to the outcome
#' @param z A numeric vector corresponding to the primary treatment or exposure
#' of interest. May be continuous or binary.
#' @param stratum A vector indicating the strata each observation belongs to.
#' @param num_trees integer number of trees to include in the BART portion of
#' the model
#' @param seed A random seed to initiate the algorithm with.
#' @param iter The total number of MCMC iterations to run.
#' @param thin Thinning rate for posterior samples.
#' @param warmup Number of MCMC iterations to run before storing posterior
#' samples.
#' @param sigma2_beta Initial value of the variance parameter used in the joint
#' multivariate normal proposal distribution for the confounding parameters.
#' This is updated during the warm-up period every `sigma2_beta_update_freq`
#' iterations.
#' @param sigma2_beta_update_freq Frequency with which to update `sigma2_beta`
#' during the warm-up period.
#' @param beta_acc_prob Target acceptance probability for the joint
#' Metropolis-Hastings proposal for the confounding parameters. `sigma2_beta`
#' will be updated throughout the MCMC with the goal of achieving this rate.
#' @param n_min Minimum number of cases for proposed terminal nodes. Moves that
#' would result in fewer than `n_min` cases in a leaf node are not considered.
#' @param moves Possible moves for the reversible-jump algorithm (DO NOT TOUCH).
#' @param move_probs Prior probabilities corresponding to the moves in `moves`.
#' @param alpha_rho,beta_rho Parameters for the tree splitting prior as defined
#' in Chipman et al. (2010).
#' @param alpha_sigma,beta_sigma Shape and scale parameters for the inverse
#' gamma prior on the terminal node variance.
#'
#' @returns `clbart` returns a list of elements of the posterior distribution
#' generated during the fitting process. The number of posterior samples is
#' K = `(iter - warmup) / thin`.
#'  \item{beta}{A K x `ncol(x)` matrix.}
#'  \item{sigma2_beta}{A length K vector (typically one value, since sigma2_beta
#'   is only updated during the warmup period).}
#'  \item{beta_acc_rate}{A length K vector of the cumulative M-H acceptance rate
#'  for the entire `beta` vector.}
#'  \item{forests}{A length K list of lists. Each internal list represents an
#'  individual tree within the ensemble, represented by a list of nodes within
#'  within that tree for the given iteration. Each node is represented as a list
#'  containing: `leaf` (indicates whether the node is a leaf node), `split_rule`
#'  (indicating the variable and cut-point defining the node), and `mu`
#'  (the prediction for the node).}
#'  \item{tree_acc_rate}{A length K vector of the proportion of tree proposals
#'  accepted out of `num_trees` within each iteration.}
#'  \item{sigma2_mu}{A length K vector of the variance parameter placed on the
#'  Normal leaf node prior.}
#'  \item{split_props}{A K x `ncol(w)` matrix containing the distribution of
#'  all splits in the ensemble across the variables in `w` at each iteration.}
#'  \item{tree_props}{A K x `ncol(w)` matrix containing the proportion of trees
#'  in the ensemble splitting on each variable in `w` at each iteration.}
#'  \item{avg_tree_depth}{A length K vector of the average depth of all trees in
#'  the ensemble at each iteration.}
#'  \item{avg_num_nodes}{A length K vector of the average number of nodes in all
#'  trees in the ensemble at each iteration.}
#'  \item{avg_num_leaves}{A length K vector of the average number of leaves in
#'  all trees in the ensemble at each iteration.}
#'  \item{time}{A length K vector of the time taken to complete each iteration.}
#'  \item{logLik}{A length K vector of the log-likelihood computed at each iteration.}
#'  \item{pWAIC, lppd, WAIC}{The WAIC penalty, log posterior predictive density,
#'  and WAIC (Widely Applicable Information Criteria) computed as:
#'  `WAIC = -2*lppd +2*pWAIC`}
#'  \item{...}{(the list also includes the parameters provided when fitting
#'  the function.)}
#'
#'
#' @examples
#' 'INSERT EXAMPLE HERE'
#'
#' @references
#' Chipman, H. A., George, E. I., and McCulloch, R. E. (2010). BART: Bayesian additive
#' regression trees. The Annals of Applied Statistics, 4(1).
#' @export
clbart <- function(w, x = NULL, y, z, stratum,
                    num_trees = 5,
                    seed = 2187,
                    iter = 1e4, thin = 1, warmup = floor(iter / 2),
                    sigma2_beta = 1,
                    sigma2_beta_update_freq = 100,
                    beta_acc_prob = 0.35,
                    n_min = 5,
                    moves = c('grow','prune','change'),
                    move_probs = c(0.3, 0.3, 0.4),
                    alpha_rho = .95, beta_rho = 2,
                    alpha_sigma = 0.0001, beta_sigma = 0.0001){

  set.seed(seed)

  # Drop columns of w that are only have one value
  w_drop <- setdiff(colnames(w), colnames(Filter(\(x) stats::var(x) > 0, w)))
  w_keep <- colnames(w)[which(!(colnames(w) %in% w_drop))]
  if(length(w_drop) > 0) print(paste(c('Dropping the following variables:', w_drop), collapse = ' '))
  w <- w[w_keep]

  # Number of strata and confounders
  stratum <- as.numeric(as.factor(stratum))
  firsts <- match(unique(stratum), stratum)
  n <- length(unique(stratum))

  # Parameters needed to accelerate computation of log-likelihood on the entire dataset
  windows <- table(stratum)
  max_win <- max(windows)
  na_locs <- (1:(max_win * n)) %in% (max_win * unique(stratum)[which(windows != max_win)]) |>
    matrix(nrow = max_win, ncol = n)

  # Set up storage for MCMC results
  K <- (iter - warmup) / thin
  k_keep <- seq(warmup + 1, iter, by = thin)

  forest_store      <- list()
  sigma2_mu_store   <- numeric(K)
  tree_acc_rate     <- numeric(K)
  split_props_store <- matrix(0, nrow = K, ncol = ncol(w), dimnames = list(NULL, colnames(w)))
  tree_props_store  <- matrix(0, nrow = K, ncol = ncol(w), dimnames = list(NULL, colnames(w)))
  avg_tree_depth    <- numeric(K)
  avg_num_nodes     <- numeric(K)
  avg_num_leaves    <- numeric(K)
  logLik_store      <- numeric(K)
  time_store        <- numeric(K)

  # Allocate running storage for WAIC calculation
  ll <- ll2 <- ell <- numeric(n)

  # Initialize beta and mu
  if(!is.null(x)){
    x <- as.matrix(x)
    p <- ncol(x)
    beta_store        <- matrix(0, ncol = p, nrow = K, dimnames = list(NULL, paste0('beta', 1:p)))
    sigma2_beta_store <- numeric(K)
    beta_acc_rate     <- numeric(K)

    m0 <- clr(y, cbind(x, z), stratum) # clogit(y ~ x + z + strata(stratum))
    beta <- m0$beta[1:p] # stats::coef(m0)[1:p]
    mu_start <- m0$beta[p+1] # stats::coef(m0)[p+1]
    beta_cov <- m0$vcov[1:p, 1:p] # stats::vcov(m0)[1:p, 1:p]
    sigma2_mu <- m0$beta[p+1] * 1e10 # stats::coef(m0)[p+1] * 1e10
    xbeta <- x %*% beta
  }
  else{
    beta_store        <- NULL
    sigma2_beta_store <- NULL
    beta_acc_rate     <- NULL
    sigma2_beta_update_freq <- NULL
    sigma2_beta <- NULL
    beta_acc_prob <- NULL

    m0 <- clr(y, matrix(z, ncol = 1), stratum) # clogit(y ~ z + strata(stratum))
    mu_start <- m0$beta
    sigma2_mu <- m0$beta * 1e10
    xbeta <- numeric(length(y))
  }

  # Initialize BART Structure
  lambda <- numeric(length(y))

  forest <- list()
  for(j in 1:num_trees){
    forest[[j]] <- plant(w, y, z, sc1 = xbeta, stratum, lambda, sigma2_mu, m_start = mu_start)
    lambda <- lambda + predict_tree(forest[[j]])
  }

  # Track beta acceptance rate (used to tune proposal distribution)
  n_acc_beta <- 0
  n_acc_beta_update <- 0

  # Begin MCMC
  cat('Beginning MCMC... \n')
  pb <- progress::progress_bar$new(
    format = "[:bar] Iteration :current/:total. Total time elapsed: :elapsedfull",
    total = iter, clear = FALSE, width = 100)
  mcmc.start <- Sys.time()
  for(k in 1:iter){

    start <- Sys.time()

    lambda <- predict_forest(forest)

    if(!is.null(x)){
      # Step 1: Update beta
      # Betas (proposal and current)
      beta_curr <- as.numeric(beta)
      beta_prop <- as.numeric(MASS::mvrnorm(1, beta, sigma2_beta * beta_cov))

      # Log-likelihoods (proposal and current)
      ll_curr <- comp_loglik(y = y, sc = x %*% beta_curr + lambda, stratum = stratum, max_win = max_win, na_locs = na_locs)
      ll_prop <- comp_loglik(y = y, sc = x %*% beta_prop + lambda, stratum = stratum, max_win = max_win, na_locs = na_locs)

      # Log-priors (proposal and current)
      lp_curr <- mvtnorm::dmvnorm(beta_curr, sigma = 1e10 * diag(p), log = TRUE)
      lp_prop <- mvtnorm::dmvnorm(beta_prop, sigma = 1e10 * diag(p), log = TRUE)

      # Acceptance ratio
      r <- min(0, ll_prop + lp_prop - ll_curr - lp_curr)

      # Update beta vector (or not)
      if(log(stats::runif(1)) < r) {
        beta <- beta_prop
        xbeta <- x %*% beta
        n_acc_beta <- n_acc_beta + 1
        n_acc_beta_update <- n_acc_beta_update + 1
      }

      # Check acceptance rate to see if need to update sigma2_beta
      if((k %% sigma2_beta_update_freq == 0) && (k <= warmup)){
        curr_rate_beta <- n_acc_beta_update / sigma2_beta_update_freq
        if(curr_rate_beta < 0.001) sigma2_beta <- sigma2_beta * 0.1
        else if(curr_rate_beta < 0.05) sigma2_beta <- sigma2_beta * 0.5
        else if(curr_rate_beta < 0.20) sigma2_beta <- sigma2_beta * 0.9
        else if(curr_rate_beta > 0.50) sigma2_beta <- sigma2_beta * 1.1
        else if(curr_rate_beta > 0.75) sigma2_beta <- sigma2_beta * 2.0
        else if(curr_rate_beta > 0.95) sigma2_beta <- sigma2_beta * 10.0
        n_acc_beta_update <- 0
      }
    }

    # Step 2: Update Tree Structure
    # Track number of trees being updated within each iteration
    n_acc_tree <- 0

    # Iterate through trees
    for(j in 1:num_trees){

      curr_tree <- forest[[j]]

      # Subtract the current tree's contribution
      lambda <- lambda - predict_tree(curr_tree)

      if(curr_tree$num_nodes == 1){
        move <- 'grow'
        p_grow <- 1
        p_prune <- 0
        p_change <- 0
      }
      else{
        p_grow <- ifelse(is_growable(curr_tree), move_probs[1], 0)
        p_prune <- move_probs[2]
        p_change <- ifelse(is_changeable(curr_tree), move_probs[3], 0)

        p_all <- sum(p_grow, p_prune, p_change)
        p_grow <- p_grow / p_all
        p_prune <- p_prune / p_all
        p_change <- p_change / p_all

        move <- sample(moves, 1, prob = c(p_grow, p_prune, p_change))
      }

      if(move == 'grow'){

        # Propose move
        prop_tree <- grow(curr_tree, w, y, z, stratum, sc1 = xbeta, lambda, sigma2_mu, n_min)
        p_grow2 <- ifelse(is_growable(prop_tree), move_probs[1], 0)
        p_prune2 <- move_probs[2]
        p_change2 <- ifelse(is_changeable(prop_tree), move_probs[3], 0)

        p_all2 <- sum(p_grow2, p_prune2, p_change2)
        p_grow2 <- p_grow2 / p_all2
        p_prune2 <- p_prune2 / p_all2
        p_change2 <- p_change2 / p_all2

        # Current and proposed nodes
        l_children_idx <- setdiff(prop_tree$node_idx, curr_tree$node_idx)
        l_idx <- get_parent(l_children_idx[1])
        l <- curr_tree$nodes[[l_idx]]
        l_children <- prop_tree$nodes[l_children_idx]

        # Current and proposed log likelihoods
        curr_logLik <- l$logLik
        prop_logLik <- sum(unlist(lapply(l_children, '[[', 'logLik')))

        # Proposal distribution priors
        g_prune <- stats::dnorm(l$mu, l$m, l$v, log = TRUE)
        g_grow <- sum(unlist(lapply(l_children, \(l) stats::dnorm(l$mu, l$m, l$v, log = TRUE))))

        # Split priors
        curr_rho <- p_split(d = l$depth, alpha = alpha_rho, beta = beta_rho)
        prop_rho <- p_split(d = l$depth + 1, alpha = alpha_rho, beta = beta_rho)

        # M-H Ratio
        r1 <- log(curr_rho) + 2 * log(1 - prop_rho) - log(1 - curr_rho)
        r2 <- prop_logLik - curr_logLik
        r3 <- log(p_prune2) - log(prop_tree$num_nogs) - log(p_grow) + log(curr_tree$num_leaves)
        r4 <- g_prune - g_grow
        r <- min(sum(r1, r2, r3, r4), 0)

        nodes_to_update <- l_children_idx
      }
      else if(move == 'prune'){

        # Propose move
        prop_tree <- prune(curr_tree, w, y, z, stratum, sc1 = xbeta, lambda, sigma2_mu, n_min)
        p_grow2 <- move_probs[1]
        p_prune2 <- ifelse(prop_tree$num_nodes != 1, move_probs[2], 0)
        p_change2 <- ifelse(is_changeable(prop_tree), move_probs[3], 0)

        p_all2 <- sum(p_grow2, p_change2, p_prune2)
        p_grow2 <- p_grow2 / p_all2
        p_prune2 <- p_prune2 / p_all2
        p_change2 <- p_change2 / p_all2

        # Current and proposed nodes
        b_children_idx <- setdiff(curr_tree$node_idx, prop_tree$node_idx)
        b_idx <- get_parent(b_children_idx[1])
        b_children <- curr_tree$nodes[b_children_idx]
        b <- prop_tree$nodes[[b_idx]]

        # Current and proposed log likelihoods
        curr_logLik <- sum(unlist(lapply(b_children, '[[', 'logLik')))
        prop_logLik <- b$logLik

        # Proposal distribution priors
        g_grow <- sum(unlist(lapply(b_children, \(l) stats::dnorm(l$mu, l$m, l$v, log = TRUE))))
        g_prune <- stats::dnorm(b$mu, b$m, b$v, log = TRUE)

        # Split priors
        curr_rho <- p_split(d = b$depth + 1, alpha = alpha_rho, beta = beta_rho)
        prop_rho <- p_split(d = b$depth, alpha = alpha_rho, beta = beta_rho)

        # M-H Ratio
        r1 <- log(1 - prop_rho) - log(prop_rho) - 2 * log(1 - curr_rho)
        r2 <- prop_logLik - curr_logLik
        r3 <- log(p_grow2) - log(prop_tree$num_leaves) - log(p_prune) + log(curr_tree$num_nogs)
        r4 <- g_grow - g_prune
        r <- min(0, sum(r1, r2, r3, r4))

        nodes_to_update <- b_idx
      }
      else if(move == 'change'){

        # Propose move
        prop_tree <- change(curr_tree, w, y, z, stratum, sc1 = xbeta, lambda, sigma2_mu, n_min)

        # Current and proposed nodes
        b_children_idx <- which(mapply(\(x, y) x != y,
                                       x = lapply(curr_tree$nodes, '[[', 'mu'),
                                       y = lapply(prop_tree$nodes, '[[', 'mu')) == 1)
        curr_b_children <- curr_tree$nodes[b_children_idx]
        prop_b_children <- prop_tree$nodes[b_children_idx]

        # Current and proposed log likelihoods
        curr_logLik <- sum(unlist(lapply(curr_b_children, '[[', 'logLik')))
        prop_logLik <- sum(unlist(lapply(prop_b_children, '[[', 'logLik')))

        # Compute the log-likelihood of the proposal
        curr_g_change <- sum(unlist(lapply(curr_b_children, \(l) stats::dnorm(l$mu, l$m, l$v, log = TRUE))))
        prop_g_change <- sum(unlist(lapply(prop_b_children, \(l) stats::dnorm(l$mu, l$m, l$v, log = TRUE))))

        # Calculate the M-H acceptance ratio on the log scale
        r <- min(prop_logLik - curr_logLik + curr_g_change - prop_g_change, 0)

        nodes_to_update <- b_children_idx
      }

      if(log(stats::runif(1)) < r){
        forest[[j]] <- prop_tree
        n_acc_tree <- n_acc_tree + 1
      }

      forest[[j]] <- update_mu(forest[[j]], y, z, stratum, xbeta, lambda, sigma2_mu)

      lambda <- lambda + predict_tree(forest[[j]])
    }

    sigma2_mu <- update_sigma2_mu(forest, alpha_sigma, beta_sigma)

    # Store Results
    if(k %in% k_keep){

      Kk <- which(k_keep == k) # ID posterior sample index

      # Fixed effects and related parameters
      if(!is.null(x)){
        beta_store[Kk,]         <- beta
        sigma2_beta_store[Kk]   <- sigma2_beta
        beta_acc_rate[Kk]       <- n_acc_beta / k
      }

      # Tree leaf node parameters
      sigma2_mu_store[Kk]     <- sigma2_mu

      # Tree structure related parameters
      forest_store[[Kk]]      <- simplify_forest(forest)
      split_props_store[Kk,]  <- get_split_props(forest, w)
      tree_props_store[Kk,]   <- get_tree_props(forest, w)
      tree_acc_rate[Kk]       <- n_acc_tree / num_trees
      avg_num_leaves[Kk]      <- get_avg_num_leaves(forest)
      avg_num_nodes[Kk]       <- get_avg_num_nodes(forest)
      avg_tree_depth[Kk]      <- get_avg_tree_depth(forest)

      # WAIC
      ll_curr <- comp_loglik(y = y, sc = xbeta + (lambda * z), stratum = stratum, max_win = max_win, na_locs = na_locs, sum = FALSE)
      ll <- ll + ll_curr
      ll2 <- ll2 + ll_curr^2
      ell <- ell + exp(ll_curr)

      # Additional info
      logLik_store[Kk]       <- sum(ll_curr)
      time_store[Kk]         <- Sys.time() - start
    }

    pb$tick()
  }

  # Return results
  pWAIC <- sum(ll2 / K - (ll / K)^2)
  lppd <- sum(log(ell / K))
  post <- list(beta = beta_store,
               sigma2_beta = sigma2_beta_store,
               beta_acc_rate = beta_acc_rate,
               forests = forest_store,
               tree_acc_rate = tree_acc_rate,
               sigma2_mu = sigma2_mu_store,
               split_props = split_props_store,
               tree_props = tree_props_store,
               avg_tree_depth = avg_tree_depth,
               avg_num_nodes = avg_num_nodes,
               avg_num_leaves = avg_num_leaves,
               time = time_store,
               logLik = logLik_store,
               pWAIC = pWAIC,
               lppd = lppd,
               WAIC = -2*lppd +2*pWAIC,
               num_trees = num_trees,
               seed = seed, iter = iter, thin = thin, warmup = warmup,
               sigma2_beta_update_freq = sigma2_beta_update_freq,
               beta_acc_prob = beta_acc_prob,
               n_min = n_min,
               moves = moves, move_probs = move_probs,
               alpha_rho = alpha_rho, beta_rho = beta_rho)

  cat('Finished! \n')

  return(post)
}
