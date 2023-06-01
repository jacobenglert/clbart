
plant <- function(w, y, z, strata, sc1 = numeric(length(y)), lambda = numeric(length(y)), sigma2_mu, m_start = 0){

  # Instantiate tree
  nodes <- list()

  # Instantiate root node (depth 0, is a leaf, not a grandparent)
  data_idx <- rep(TRUE, length(y))
  data_idx_first <- data_idx & (1:length(y) %in% match(unique(strata), strata))

  root <- list(depth = 0, leaf = 1, nog = 0,
               data_idx = data_idx,
               data_idx_first = data_idx_first,
               split_var = NULL, split_rule = NULL)
  root[c('n','max_win','na_locs')] <- get_max_win_and_na_locs(strata)

  # Compute m and v for the root node
  root[c('m','v')] <- comp_mv(m_start, y, z, strata, sc1, lambda, sigma2_mu,
                              max_win = root$max_win, na_locs = root$na_locs)

  # Sample mu for the root node and update the log-likelihood
  root['mu'] <- stats::rnorm(1, root$m, root$v)
  root['logLik'] <- comp_loglik(y, sc = sc1 + z * (lambda + root$mu), strata,
                                max_win = root$max_win, na_locs = root$na_locs)

  nodes[[1]] <- root

  return(list(nodes = nodes, node_idx = 1, num_nodes = 1, num_leaves = 1,
              num_nogs = 0, max_depth = 1,
              valid_grow_nodes = 1, valid_grow_vars = list(colnames(w)),
              valid_change_nodes = numeric(0), valid_change_vars = list()))
}

grow <- function(tree, w, y, z, strata, sc1 = numeric(length(y)), lambda = numeric(length(y)), sigma2_mu, n_min = 1){

  stopifnot("Cannot grow additional leaves due to lack of available splits" = length(tree$valid_grow_nodes) != 0)

  l_idx <- resample(tree$valid_grow_nodes, 1)
  l <- tree$nodes[[l_idx]]
  var_s <- sample(tree$valid_grow_vars[[which(tree$valid_grow_nodes == l_idx)]], 1)
  var_c <- resample(get_valid_split_rules(x = w[l$data_idx_first, var_s, drop = TRUE], n_min = n_min), 1)

  id_left <- (w[, var_s, drop = TRUE] <= var_c) & l$data_idx
  id_right <- (w[, var_s, drop = TRUE] > var_c) & l$data_idx
  child_data_idx <- list(id_left, id_right)
  splits <- list(paste(var_s, '<=', var_c), paste(var_s, '>', var_c))

  tree$nodes[[l_idx]][c('leaf','nog')] <- c(0, 1)

  if(l_idx != 1) tree$nodes[[get_parent(l_idx)]]$nog <- 0

  l_children <- get_children(l_idx)
  for(i in 1:2){

    idx <- l_children[i]

    data_idx_first <- child_data_idx[[i]] & (1:length(strata) %in% match(unique(strata), strata))

    child <- list(depth = l$depth + 1, leaf = 1, nog = 0,
                  data_idx = child_data_idx[[i]],
                  data_idx_first = data_idx_first)
    child[c('n','max_win','na_locs')] <- get_max_win_and_na_locs(strata[child$data_idx])
    child[c('m','v')] <- comp_mv(m_start = l$m,
                                 y = y[child$data_idx],
                                 z = z[child$data_idx],
                                 strata = strata[child$data_idx],
                                 sc1 = sc1[child$data_idx],
                                 lambda = lambda[child$data_idx],
                                 sigma2_mu = sigma2_mu,
                                 max_win = child$max_win,
                                 na_locs = child$na_locs)

    child['mu'] <- stats::rnorm(1, child$m, child$v)
    child['logLik'] <- comp_loglik(y = y[child$data_idx],
                                   sc = (sc1 + z * (lambda + child$mu))[child$data_idx],
                                   strata = strata[child$data_idx],
                                   max_win = child$max_win,
                                   na_locs = child$na_locs)

    child$split_var <- var_s
    child$split_rule <- c(l$split_rule, splits[[i]])

    tree$nodes[[idx]] <- child
  }

  tree <- fixTreeParms(tree)
  tree[c('valid_grow_vars', 'valid_grow_nodes')] <- get_grow_moves(tree, w, n_min)
  tree[c('valid_change_vars', 'valid_change_nodes')] <- get_change_moves(tree, w, n_min)

  return(tree)
}

prune <- function(tree, w, y, z, strata, sc1 = numeric(length(y)), lambda = numeric(length(y)), sigma2_mu, n_min = 1){

  stopifnot("Cannot prune a tree that only contains the root node" = tree$num_nodes != 1)

  b_idx <- resample(which(lapply(tree$nodes, \(n) as.numeric(n$nog)) == 1), 1)

  if((b_idx != 1) && (tree$nodes[[get_sibling(b_idx)]]$leaf == 1)){
    tree$nodes[[get_parent(b_idx)]]$nog <- 1
  }

  b <- tree$nodes[[b_idx]]
  b[c('leaf','nog')] <- c(1, 0)

  b[c('m','v')] <- comp_mv(m_start = sum(unlist(lapply(tree$nodes[get_children(b_idx)], \(l) l$n * l$m))) / b$n,
                           y = y[b$data_idx],
                           z = z[b$data_idx],
                           strata = strata[b$data_idx],
                           sc1 = sc1[b$data_idx],
                           lambda = lambda[b$data_idx],
                           sigma2_mu = sigma2_mu,
                           max_win = b$max_win,
                           na_locs = b$na_locs)

  b['mu'] <- stats::rnorm(1, b$m, b$v)
  b['logLik'] <- comp_loglik(y = y[b$data_idx],
                             sc = (sc1 + z * (lambda + b$mu))[b$data_idx],
                             strata = strata[b$data_idx],
                             max_win = b$max_win,
                             na_locs = b$na_locs)

  tree$nodes[[b_idx]] <- b

  tree$nodes[get_children(b_idx)] <- list(NULL)

  tree <- fixTreeParms(tree)
  tree[c('valid_grow_vars', 'valid_grow_nodes')] <- get_grow_moves(tree, w, n_min)
  tree[c('valid_change_vars', 'valid_change_nodes')] <- get_change_moves(tree, w, n_min)

  return(tree)
}

change <- function(tree, w, y, z, strata, sc1 = numeric(length(y)), lambda = numeric(length(y)), sigma2_mu, n_min = 1){

  stopifnot("Cannot change a tree that only contains the root node" = tree$num_nodes != 1)

  stopifnot("Cannot change rule due to lack of available splits" = length(tree$valid_change_nodes) != 0)

  b_idx <- resample(tree$valid_change_nodes, 1)
  b <- tree$nodes[[b_idx]]
  var_s <- sample(tree$valid_change_vars[[which(tree$valid_change_nodes == b_idx)]], 1)
  var_c <- resample(get_valid_split_rules(x = w[b$data_idx_first, var_s, drop = TRUE], n_min = n_min), 1)

  id_left <- (w[, var_s, drop = TRUE] <= var_c) & b$data_idx
  id_right <- (w[, var_s, drop = TRUE] > var_c) & b$data_idx
  child_data_idx <- list(id_left, id_right)
  splits <- list(paste(var_s, '<=', var_c), paste(var_s, '>', var_c))

  b_children <- get_children(b_idx)
  for(i in 1:2){

    idx <- b_children[i]
    data_idx_first <- child_data_idx[[i]] & (1:length(strata) %in% match(unique(strata), strata))

    child <- list(depth = b$depth + 1, leaf = 1, nog = 0,
                  data_idx = child_data_idx[[i]],
                  data_idx_first = data_idx_first)
    child[c('n','max_win','na_locs')] <- get_max_win_and_na_locs(strata[child$data_idx])
    child[c('m','v')] <- comp_mv(m_start = b$m,
                                 y = y[child$data_idx],
                                 z = z[child$data_idx],
                                 strata = strata[child$data_idx],
                                 sc1 = sc1[child$data_idx],
                                 lambda = lambda[child$data_idx],
                                 sigma2_mu = sigma2_mu,
                                 max_win = child$max_win,
                                 na_locs = child$na_locs)

    child['mu'] <- stats::rnorm(1, child$m, child$v)
    child['logLik'] <- comp_loglik(y = y[child$data_idx],
                                   sc = (sc1 + z * (lambda + child$mu))[child$data_idx],
                                   strata = strata[child$data_idx],
                                   max_win = child$max_win,
                                   na_locs = child$na_locs)

    child$split_var <- var_s
    child$split_rule <- c(b$split_rule, splits[[i]])

    tree$nodes[[idx]] <- child
  }

  tree <- fixTreeParms(tree)
  tree[c('valid_grow_vars', 'valid_grow_nodes')] <- get_grow_moves(tree, w, n_min)
  tree[c('valid_change_vars', 'valid_change_nodes')] <- get_change_moves(tree, w, n_min)

  return(tree)
}
