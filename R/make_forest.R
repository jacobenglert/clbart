# make_forest <- function(w, y, z, s, offset, hypers) {
# 
#   # root_fit <- get_node_map(z, y, s, offset, 1:nrow(w), hypers)
#   # m <- as.numeric(root_fit$m)
#   # v <- as.numeric(root_fit$v)
# 
#   forest <- data.frame(tree = 1:hypers$num_tree)
#   forest$node <- 1
#   forest$n <- nrow(w)
#   forest$var <- -1
#   forest$nog <- 0
#   forest$m <- NA # m
#   forest$v <- NA # v
#   forest$value <- 0 # rnorm(hypers$num_tree, forest$m, forest$v)
#   for(t in 1:hypers$num_tree) forest$loglik[t] <- as.numeric(clr_loglik_cpp(beta = forest$value[t],
#                                              x = z,
#                                              y = y,
#                                              stratum = s,
#                                              offset = offset))
#   forest$goodvars <- list(get_goodvars(w))
#   return(forest)
# }


make_forest <- function(w, y, z, s, offset, hypers) {

  # Create shell
  forest <- data.frame(tree = 1:hypers$num_tree)
  forest$node <- 1
  forest$n    <- nrow(w)
  forest$var  <- -1
  forest$nog  <- 0

  # # Get initial fit
  # init.fit <- clr_cpp(cbind(x, z), y, strata)
  # beta <- init.fit$beta[length(init.fit$beta)]
  # lambda <- rep(beta, nrow(w))
  
  # Initialize values with generalized backfitting
  for (t in 1:hypers$num_tree) {

    root_fit <- get_node_map(z, y, s, offset, 1:nrow(w), hypers)
    m <- root_fit$m
    v <- root_fit$v

    forest$m[t] <- m
    forest$v[t] <- v
    forest$value[t]   <- stats::rnorm(1, m, v)
    forest$loglik[t]  <- as.numeric(clr_loglik_cpp(beta = forest$value[t],
                                        x = z,
                                        y = y,
                                        stratum = s,
                                        offset = offset))
    
    offset <- offset + forest$value[t] * z

  }

  forest$goodvars <- list(get_goodvars(w))

  return(forest)
}
