# Grow proposal
grow <- function (tree, w, y, z, s, offset, hypers) {

  # Store current tree
  curr_tree <- tree

  # Retrieve leaf nodes from current tree
  leaves <- curr_tree[curr_tree$var == -1,]

  # Exclude leaf nodes that cannot be grown any further
  can_grow <- sapply(leaves$goodvars, \(l) any(unlist(l)))
  leaves <- leaves[can_grow,]

  if (nrow(leaves) == 0) return(tree)

  # Sample a leaf node
  l <- leaves[sample(nrow(leaves), size = 1),]

  # Identify observations mapped to this node
  l.ind <- get_mapped_observations(curr_tree, l$node, w)


  # Update the current tree (selected leaf only) --------------------------
  # Necessary for the acceptance ratio to be calculated correctly
  l.fit <- get_node_map(z, y, s, offset, l.ind, hypers, start_value = 0) # na.omit(c(l$m, 0))[1])
  l$m <- l.fit$m
  l$v <- l.fit$v
  l$loglik <- as.numeric(clr_loglik_cpp(beta = l$value, # uses same value from before
                             x = z[l.ind,,drop = FALSE],
                             y = y[l.ind],
                             stratum = match(s[l.ind], unique(s[l.ind])),
                             offset = offset[l.ind]))
  curr_tree[curr_tree$node == l$node,] <- l


  # Grow the current tree at the selected leaf node -----------------------

  # Copy current tree
  prop_tree <- curr_tree

  # Sample a splitting rule
  split.probs <- hypers$s * unlist(l$goodvars)
  l$var <- sample(1:ncol(w), size = 1, prob = split.probs / sum(split.probs))
  w.sub <- w[l.ind,l$var]
  l$value <- sample(w.sub[w.sub < max(w.sub)], size = 1) # sample(w.sub[floor(runif(1) * l$n)])

  # Update existing leaf
  l$nog <- 1

  # Update parent of selected leaf to not a NOG (may have already been)
  prop_tree$nog[prop_tree$node == get_parent(l$node)] <- 0

  # Create new leaves
  lL <- lR <- l
  lL$node <- l$node * 2
  lL$var <- -1
  lL$nog <- 0
  lL.ind <- l.ind & w[,l$var] <= l$value
  lL$n <- sum(lL.ind)

  lR$node <- l$node * 2 + 1
  lR$var <- -1
  lR$nog <- 0
  lR.ind <- l.ind & (!lL.ind) # l.ind & w[,w.var] > l$value
  lR$n <- sum(lR.ind)

  if (min(lL$n, lR$n) < 20) return (curr_tree)

  # Sample node values using Laplace approximation
  lL.fit <- get_node_map(z, y, s, offset, lL.ind, hypers, start_value = l$m)
  lL$m <- lL.fit$m
  lL$v <- lL.fit$v
  if (is.nan(lL$v) | is.infinite(lL$v)) return(curr_tree)
  lL$value <- stats::rnorm(1, lL$m, lL$v)
  lL$loglik <- as.numeric(clr_loglik_cpp(beta = lL$value,
                              x = z[lL.ind,,drop = FALSE],
                              y = y[lL.ind],
                              stratum = match(s[lL.ind], unique(s[lL.ind])),
                              offset = offset[lL.ind]))
  lL$goodvars <- list(get_goodvars(w, lL.ind))

  lR.fit <- get_node_map(z, y, s, offset, lR.ind, hypers, start_value = l$m)
  lR$m <- lR.fit$m
  lR$v <- lR.fit$v
  if (is.nan(lR$v) | is.infinite(lR$v)) return(curr_tree)
  lR$value <- stats::rnorm(1, lR$m, lR$v)
  lR$loglik <- as.numeric(clr_loglik_cpp(beta = lR$value,
                              x = z[lR.ind,,drop = FALSE],
                              y = y[lR.ind],
                              stratum = match(s[lR.ind], unique(s[lR.ind])),
                              offset = offset[lR.ind]))
  lR$goodvars <- list(get_goodvars(w, lR.ind))

  # Update proposed tree
  prop_tree <- rbind(prop_tree[prop_tree$node != l$node,], l, lL, lR)
  prop_tree <- prop_tree[order(prop_tree$node),]


  # Accept/Reject Proposal ------------------------------------------------

  # Compute M-H acceptance ratio on log scale
  # r <- compute_grow_ratio(old_tree = curr_tree, new_tree = prop_tree, hypers)

  # Retrieve current parameter for leaf node
  l_value <- curr_tree$value[curr_tree$node == l$node]

  # Prior ratio
  r_prior <- log(node_depth_prior(l$node, alpha = hypers$base, beta = hypers$power)) +
    2 * log(1 - node_depth_prior(lL$node, alpha = hypers$base, beta = hypers$power)) -
    log(1 - node_depth_prior(l$node, alpha = hypers$base, beta = hypers$power)) +
    stats::dnorm(lL$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) +
    stats::dnorm(lR$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) -
    stats::dnorm(l_value, hypers$mu_mu, hypers$sigma_mu, log = TRUE)

  # Likelihood ratio
  r_loglik <- lL$loglik + lR$loglik - l$loglik

  # Proposal ratio
  can_grow <- sapply(curr_tree$goodvars, \(l) any(unlist(l)))
  r_prop <- log(p_moves(prop_tree)['p_prune']) -
    log(sum(prop_tree$nog)) -
    log(p_moves(curr_tree)['p_grow']) +
    log(sum((curr_tree$var == -1) & can_grow)) +
    stats::dnorm(l_value, l$m, l$v, log = TRUE) -
    stats::dnorm(lL$value, lL$m, lL$v, log = TRUE) -
    stats::dnorm(lR$value, lR$m, lR$v, log = TRUE)

  r <- r_prior + r_loglik + r_prop

  if (r > log(stats::runif(1))) {
    return(prop_tree)
  } else {
    return(curr_tree)
  }

}


# Prune proposal
prune <- function(tree, w, y, z, s, offset, hypers){

  # Store current tree
  curr_tree <- tree

  # Return original tree if only consists of the root node
  if (nrow(tree) == 1) return(tree)

  # Sample a NOG node
  nogs <- curr_tree[curr_tree$nog == 1,]
  b <- nogs[sample(nrow(nogs), 1),]

  # Identify observations mapped to this node
  b.ind <- get_mapped_observations(curr_tree, b$node, w)

  # Get leaves beneath this node
  bL <- curr_tree[curr_tree$node == get_children(b$node)[1],]
  bR <- curr_tree[curr_tree$node == get_children(b$node)[2],]

  # Update the current tree (selected leaves only) ------------------------
  # Necessary for the acceptance ratio to be calculated correctly
  bL.ind <- b.ind & (w[,b$var] <= b$value)
  bL.fit <- get_node_map(z, y, s, offset, bL.ind, hypers, start_value = bL$m)
  bL$m <- bL.fit$m
  bL$v <- bL.fit$v
  bL$loglik <- as.numeric(clr_loglik_cpp(beta = bL$value,
                              x = z[bL.ind,,drop = FALSE],
                              y = y[bL.ind],
                              stratum = match(s[bL.ind], unique(s[bL.ind])),
                              offset = offset[bL.ind]))

  bR.ind <- b.ind & !bL.ind
  bR.fit <- get_node_map(z, y, s, offset, bR.ind, hypers, start_value = bR$m)
  bR$m <- bR.fit$m
  bR$v <- bR.fit$v
  bR$loglik <- as.numeric(clr_loglik_cpp(beta = bR$value,
                              x = z[bR.ind,,drop = FALSE],
                              y = y[bR.ind],
                              stratum = match(s[bR.ind], unique(s[bR.ind])),
                              offset = offset[bR.ind]))

  curr_tree[curr_tree$node == bL$node,] <- bL
  curr_tree[curr_tree$node == bR$node,] <- bR


  # Prune tree back to the selected NOG node ------------------------------

  # Copy current tree
  prop_tree <- curr_tree

  # Remove the leaves beneath the selected NOG node
  prop_tree <- prop_tree[!(prop_tree$node %in% get_children(b$node)),]

  # Convert the selected NOG node to a leaf node
  b$var <- -1
  b$nog <- 0

  # Convert parent of b to a NOG node if b's sibling is already a leaf
  if (b$node != 1) {
    if (prop_tree$var[prop_tree$node == get_sibling(b$node)] == -1) {
      prop_tree$nog[prop_tree$node == get_parent(b$node)] <- 1
    }
  }

  # Sample node value using Laplace approximation
  b.fit <- get_node_map(z, y, s, offset, b.ind, hypers, start_value = (bL$n * bL$m + bR$n * bR$m) / (bL$n + bR$n))
  b$m <- b.fit$m
  b$v <- b.fit$v
  if (is.nan(b$v) | is.infinite(b$v)) return(curr_tree)
  b$value <- stats::rnorm(1, b$m, b$v)
  b$loglik <- as.numeric(clr_loglik_cpp(beta = b$value,
                             x = z[b.ind,,drop = FALSE],
                             y = y[b.ind],
                             stratum = match(s[b.ind], unique(s[b.ind])),
                             offset = offset[b.ind]))

  # Update proposed tree
  prop_tree[prop_tree$node == b$node,] <- b
  prop_tree <- prop_tree[order(prop_tree$node),]


  # Accept/Reject Proposal ------------------------------------------------

  # Compute M-H acceptance ratio on log scale
  # r <- compute_prune_ratio(old_tree = curr_tree, new_tree = prop_tree, hypers)

  # Prior ratio
  r_prior <- log(1 - node_depth_prior(b$node, alpha = hypers$base, beta = hypers$power)) -
    log(node_depth_prior(b$node, alpha = hypers$base, beta = hypers$power)) -
    2 * log(1 - node_depth_prior(bL$node, alpha = hypers$base, beta = hypers$power)) +
    stats::dnorm(b$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) -
    stats::dnorm(bL$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) -
    stats::dnorm(bR$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE)

  # Likelihood ratio
  r_loglik <- b$loglik - bL$loglik - bR$loglik

  # Proposal ratio
  can_grow <- sapply(prop_tree$goodvars, \(l) any(unlist(l)))
  r_prop <- log(p_moves(prop_tree)['p_grow']) -
    log(sum((prop_tree$var == -1) & can_grow)) -
    log(p_moves(curr_tree)['p_prune']) +
    log(sum(curr_tree$nog)) +
    stats::dnorm(bL$value, bL$m, bL$v, log = TRUE) +
    stats::dnorm(bR$value, bR$m, bR$v, log = TRUE) -
    stats::dnorm(b$value, b$m, b$v, log = TRUE)

  r <- r_prior + r_loglik + r_prop

  if (r > log(stats::runif(1))) {
    return(prop_tree)
  } else {
    return(curr_tree)
  }
}


# Change proposal
change <- function(tree, w, y, z, s, offset, hypers){

  # Store current tree
  curr_tree <- tree

  # Return original tree if only consists of the root node
  if (nrow(tree) == 1) return(tree)

  # Retrieve NOG nodes
  nogs <- curr_tree[curr_tree$nog == 1,]

  # Sample a NOG node
  b <- nogs[sample(nrow(nogs), 1),]

  # Identify observations mapped to this node
  b.ind <- get_mapped_observations(curr_tree, b$node, w)

  # Get leaves beneath this node
  bL <- curr_tree[curr_tree$node == get_children(b$node)[1],]
  bR <- curr_tree[curr_tree$node == get_children(b$node)[2],]

  # Update the current tree (selected leaves only) ------------------------
  # Necessary for the acceptance ratio to be calculated correctly
  bL.ind <- b.ind & (w[,b$var] <= b$value)
  bL.fit <- get_node_map(z, y, s, offset, bL.ind, hypers, start_value = bL$m)
  bL$m <- bL.fit$m
  bL$v <- bL.fit$v
  bL$loglik <- as.numeric(clr_loglik_cpp(beta = bL$value,
                              x = z[bL.ind,,drop = FALSE],
                              y = y[bL.ind],
                              stratum = match(s[bL.ind], unique(s[bL.ind])),
                              offset = offset[bL.ind]))

  bR.ind <- b.ind & !bL.ind
  bR.fit <- get_node_map(z, y, s, offset, bR.ind, hypers, start_value = bR$m)
  bR$m <- bR.fit$m
  bR$v <- bR.fit$v
  bR$loglik <- as.numeric(clr_loglik_cpp(beta = bR$value,
                              x = z[bR.ind,,drop = FALSE],
                              y = y[bR.ind],
                              stratum = match(s[bR.ind], unique(s[bR.ind])),
                              offset = offset[bR.ind]))

  curr_tree[curr_tree$node == bL$node,] <- bL
  curr_tree[curr_tree$node == bR$node,] <- bR


  # Change the splitting criteria at the selected NOG node ----------------

  # Copy current tree
  prop_tree <- curr_tree

  # Sample a valid splitting rule and update selected NOG node
  split.probs <- hypers$s * unlist(b$goodvars)
  b$var <- sample(1:ncol(w), size = 1, prob = split.probs / sum(split.probs))
  w.sub <- w[b.ind,b$var]
  b$value <- sample(w.sub[w.sub < max(w.sub)], size = 1)

  # Propose new leaves for selected NOG node
  bL.prop <- bL
  bR.prop <- bR

  bL.prop.ind <- b.ind & (w[,b$var] <= b$value)
  bL.prop$n <- sum(bL.prop.ind)

  bR.prop.ind <- b.ind & !bL.prop.ind
  bR.prop$n <- sum(bR.prop.ind)

  if (min(bL.prop$n, bR.prop$n) < 20) return (curr_tree)

  # Sample node values using Laplace approximation
  bL.prop.fit <- get_node_map(z, y, s, offset, bL.prop.ind, hypers, start_value = b$m)
  bL.prop$m <- bL.prop.fit$m
  bL.prop$v <- bL.prop.fit$v
  if (is.nan(bL.prop$v) | is.infinite(bL.prop$v)) return(curr_tree)
  bL.prop$value <- stats::rnorm(1, bL.prop$m, bL.prop$v)
  bL.prop$loglik <- as.numeric(clr_loglik_cpp(beta = bL.prop$value,
                                   x = z[bL.prop.ind,,drop = FALSE],
                                   y = y[bL.prop.ind],
                                   stratum = match(s[bL.prop.ind], unique(s[bL.prop.ind])),
                                   offset = offset[bL.prop.ind]))
  bL.prop$goodvars <- list(get_goodvars(w, bL.prop.ind))

  bR.prop.fit <- get_node_map(z, y, s, offset, bR.prop.ind, hypers, start_value = b$m)
  bR.prop$m <- bR.prop.fit$m
  bR.prop$v <- bR.prop.fit$v
  if (is.nan(bR.prop$v) | is.infinite(bR.prop$v)) return(curr_tree)
  bR.prop$value <- stats::rnorm(1, bR.prop$m, bR.prop$v)
  bR.prop$loglik <- as.numeric(clr_loglik_cpp(beta = bR.prop$value,
                                   x = z[bR.prop.ind,,drop = FALSE],
                                   y = y[bR.prop.ind],
                                   stratum = match(s[bR.prop.ind], unique(s[bR.prop.ind])),
                                   offset = offset[bR.prop.ind]))
  bR.prop$goodvars <- list(get_goodvars(w, bR.prop.ind))

  # Update proposed tree
  prop_tree[prop_tree$node == b$node,] <- b
  prop_tree[prop_tree$node == bL$node,] <- bL.prop
  prop_tree[prop_tree$node == bR$node,] <- bR.prop
  prop_tree <- prop_tree[order(prop_tree$node),]


  # Accept/Reject Proposal --------------------------------------------------

  # Compute M-H acceptance ratio on log scale
  # r <- compute_change_ratio(old_tree = curr_tree, new_tree = prop_tree, hypers)

  # Prior ratio
  r_prior <- stats::dnorm(bL.prop$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) +
    stats::dnorm(bR.prop$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) -
    stats::dnorm(bL$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) -
    stats::dnorm(bR$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE)

  # Likelihood ratio
  r_loglik <- bL.prop$loglik + bR.prop$loglik - bL$loglik - bR$loglik

  # Proposal ratio
  r_prop <- stats::dnorm(bL$value, bL$m, bL$v, log = TRUE) +
    stats::dnorm(bR$value, bR$m, bR$v, log = TRUE) -
    stats::dnorm(bL.prop$value, bL.prop$m, bL.prop$v, log = TRUE) -
    stats::dnorm(bR.prop$value, bR.prop$m, bR.prop$v, log = TRUE)

  r <- r_prior + r_loglik + r_prop

  if (r > log(stats::runif(1))) {
    return(prop_tree)
  } else {
    return(curr_tree)
  }
}
