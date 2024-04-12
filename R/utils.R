get_sibling <- function(x){
  if(x == 1) return(numeric(0))
  if(x %% 2 == 0) return(x + 1) else return(x-1)
}

get_children <- function(x){
  return(c(x * 2, x * 2 + 1))
}

get_parent <- function(x){
  if(x == 1) return(numeric(0))
  else return(floor(x / 2))
}

get_ancestors <- function(x){
  if(x == 1) return(numeric(0))
  else {
    a <- get_parent(x)
    while(a[length(a)] != 1) a <- c(a, get_parent(a[length(a)]))
    return(rev(a))
  }
}

get_node_depth <- function(x) return(floor(log2(x)))

node_depth_prior <- function(x, alpha = 0.95, beta = 2) return(alpha * (1 + get_node_depth(x))^(-beta))

p_grow <- function(tree, default = 1 / 3) return(ifelse(nrow(tree) == 1, 1, default))

p_prune <- function(tree, default = 1 / 3) return(ifelse(nrow(tree) == 1, 0, default))

p_moves <- function (tree, default = c(p_grow = 0.3, p_prune = 0.3, p_change = 0.4)) {

  # Check if tree can grow
  if (any(unlist(tree$goodvars[tree$var == -1]))) can_grow <- TRUE else can_grow <- FALSE

  # Check if tree can be pruned
  if (nrow(tree) > 1) can_prune <- TRUE else can_prune <- FALSE

  # If tree can be pruned it can also likely be changed (might rarely cause an issue)
  if (can_prune) can_change <- TRUE else can_change <- FALSE

  probs <- default * c(can_grow, can_prune, can_change)

  return (probs / sum(probs))
}

get_goodvars <- function (x, ind = 1:nrow(x)) {
  goodvars <- list()
  for(p in 1:ncol(x)) {
    xip <- x[ind,p]
    xip_min <- min(xip)
    xip_max <- max(xip)
    goodvars[[p]] <- xip_min != xip_max
  }
  return(goodvars)

}


log_sum_exp <- function(log_values) {
  # Find the maximum log value
  max_log_value <- max(log_values)

  # Compute the sum of exponentiated log values, adjusted by the max log value for numerical stability
  sum_exp <- sum(exp(log_values - max_log_value))

  # Return the log of this sum, re-adjusting by adding the max log value back
  return(max_log_value + log(sum_exp))
}

logs_to_s <- function(log_probs) {

  # Step 1: Subtract the maximum log probability
  max_log_prob <- max(log_probs)
  shifted_probs <- log_probs - max_log_prob

  # Step 2: Exponentiate
  probs <- exp(shifted_probs)

  # Step 3: Normalize
  probs <- probs / sum(probs)

  return(probs)

}

# rdirichlet <- function(n, alpha) {
#   k <- length(alpha)
#   dirichlet_samples <- matrix(nrow = n, ncol = k)
#
#   for (i in 1:n) {
#     gamma_samples <- rgamma(k, shape = alpha, rate = 1)
#     dirichlet_samples[i, ] <- gamma_samples / sum(gamma_samples)
#   }
#
#   return(dirichlet_samples)
# }
