# Compute the transition ratio for a GROW move on the log scale
compute_grow_ratio <- function(old_tree, new_tree, hypers) {

  # If trees are equal return 0
  if (nrow(old_tree) == nrow(new_tree)) {
    if (all(old_tree$value == new_tree$value)) return(-Inf)
  }

  # Identify old and new nodes
  shared_nodes <- intersect(old_tree$node, new_tree$node)
  new_nodes <- new_tree[!(new_tree$node %in% shared_nodes),]
  old_node <- old_tree[old_tree$node == get_parent(new_nodes$node[1]),]

  # Part 1: node depth prior ratio
  r1 <- log(node_depth_prior(old_node$node, hypers$base, hypers$power)) +
    2*(log(1 - node_depth_prior(new_nodes$node[1], hypers$base, hypers$power))) -
    log(1 - node_depth_prior(old_node$node, hypers$base, hypers$power))

  # Part 2: node value prior ratio
  r2 <- sum(stats::dnorm(new_nodes$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE)) -
    stats::dnorm(old_node$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE)

  # Part 3: likelihood ratio
  r3 <- sum(new_nodes$loglik) - old_node$loglik

  # Part 4: tree proposal prior ratio
  r4 <- log(p_moves(new_tree)['p_prune']) -
    log(sum(new_tree$nog)) -
    log(p_moves(old_tree)['p_grow']) +
    log(sum((old_tree$var == -1)[sapply(old_tree$goodvars, \(l) any(unlist(l)))]))

  # Part 5: node value proposal ratio
  r5 <- stats::dnorm(old_node$value, old_node$m, old_node$v, log = TRUE) -
    sum(stats::dnorm(new_nodes$value, new_nodes$m, new_nodes$v, log = TRUE))

  r <- min(sum(r1, r2, r3, r4, r5), 0)

  return(r)
}

# Compute the transition ratio for a PRUNE move on the log scale
compute_prune_ratio <- function(old_tree, new_tree, hypers) {

  # Identify old and new nodes
  old_nodes <- old_tree[old_tree$node %in% setdiff(old_tree$node, new_tree$node),]
  new_node <- new_tree[new_tree$node == get_parent(old_nodes$node[1]),]

  # Part 1: node depth prior ratio
  r1 <- log(1 - node_depth_prior(new_node$node, hypers$base, hypers$power)) -
    log(node_depth_prior(new_node$node, hypers$base, hypers$power)) -
    2*log(1 - node_depth_prior(old_nodes$node[1], hypers$base, hypers$power))


  # Part 2: node value prior ratio
  r2 <- stats::dnorm(new_node$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE) -
    sum(stats::dnorm(old_nodes$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE))

  # Part 3: likelihood ratio
  r3 <- new_node$loglik - sum(old_nodes$loglik)

  # Part 4: tree proposal prior ratio
  r4 <- log(p_moves(new_tree)['p_grow']) -
    log(sum((new_tree$var == -1)[sapply(new_tree$goodvars, \(l) any(unlist(l)))])) -
    log(p_moves(old_tree)['p_prune']) +
    log(sum(old_tree$nog))

  # Part 5: node value proposal ratio
  r5 <- sum(stats::dnorm(old_nodes$value, old_nodes$m, old_nodes$v, log = TRUE)) -
    stats::dnorm(new_node$value, new_node$m, new_node$v, log = TRUE)

  r <- min(sum(r1, r2, r3, r4, r5), 0)

  return(r)
}

# Compute the transition ratio for a CHANGE move on the log scale
compute_change_ratio <- function(old_tree, new_tree, hypers) {

  # If trees are equal return 0
  if (all(old_tree$value == new_tree$value)) return(-Inf)

  # Identify old and new nodes
  changed_nodes <- old_tree[old_tree$value != new_tree$value & old_tree$var == -1,]$node
  old_nodes <- old_tree[old_tree$node %in% changed_nodes,]
  new_nodes <- new_tree[new_tree$node %in% changed_nodes,]

  # Part 1: node value prior ratio
  r1 <- sum(stats::dnorm(new_nodes$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE)) -
    sum(stats::dnorm(old_nodes$value, hypers$mu_mu, hypers$sigma_mu, log = TRUE))

  # Part 2: likelihood ratio
  r2 <- sum(new_nodes$loglik) - sum(old_nodes$loglik)

  # Part 3: node value proposal ratio
  r3 <- sum(stats::dnorm(old_nodes$value, old_nodes$m, old_nodes$v, log = TRUE)) -
    sum(stats::dnorm(new_nodes$value, new_nodes$m, new_nodes$v, log = TRUE))

  r <- min(sum(r1, r2, r3), 0)

  return(r)
}

