get_max_win_and_na_locs <- function(strata){
  strata <- as.numeric(as.factor(strata))
  n <- length(unique(strata))
  windows <- table(strata)
  max_win <- max(windows)
  na_locs <- (1:(max_win * n)) %in% (max_win * unique(strata)[which(windows != max_win)]) |>
    matrix(nrow = max_win, ncol = n)
  return(list(n = n, max_win = max_win, na_locs = na_locs))
}

get_valid_split_rules <- function(x, n_sample = length(x), n_min){

  q1 <- quantile(x, n_min / n_sample, type = 1)
  q2 <- quantile(x, (n_sample - n_min + 1)/ n_sample, type = 1)

  var_range <- sort(unique(x[x >= q1 & x < q2]))

  return(var_range)
}

get_valid_split_vars <- function(x, n_sample = nrow(x), n_min){
  valid_splits <- lapply(x, function(x0){
    q1 <- quantile(x0, n_min / n_sample, type = 1)
    q2 <- quantile(x0, (n_sample - n_min + 1)/ n_sample, type = 1)
    if(n_min > n_sample / 2) FALSE
    else if(q1 != q2) TRUE
    else FALSE
  })

  return(colnames(x)[which(unlist(valid_splits))])
}

get_grow_moves <- function(tree, w, n_min){

  # ID leaves
  l_idx <- which(lapply(tree$nodes, \(n) as.numeric(n$leaf)) == 1)
  l <- tree$nodes[l_idx]

  valid_split_vars <- lapply(l, \(l) get_valid_split_vars(x = w[l$data_idx_first,], n_min = n_min))
  valid_split_nodes <- which(lapply(valid_split_vars, length) > 0)

  return(list(valid_grow_vars = valid_split_vars[valid_split_nodes],
              valid_grow_nodes = l_idx[valid_split_nodes]))
}

get_change_moves <- function(tree, w, n_min){

  # ID nogs
  b_idx <- which(lapply(tree$nodes, \(n) as.numeric(n$nog)) == 1)

  if(length(b_idx) == 0){
    return(list(valid_change_vars = list(NULL),
                valid_change_nodes = numeric(0)))
  }
  else{

    valid_split_vars <- sapply(b_idx, \(b){
      child_rule <- tree$nodes[[get_children(b)[1]]]$split_var
      valid_split_vars <- get_valid_split_vars(x = w[tree$nodes[[b]]$data_idx_first,], n_min = n_min)
      valid_split_vars <- setdiff(valid_split_vars, child_rule)
    }, simplify = FALSE)

    valid_split_nodes <- which(lapply(valid_split_vars, length) > 0)

    return(list(valid_change_vars = valid_split_vars[valid_split_nodes],
                valid_change_nodes = b_idx[valid_split_nodes]))
  }
}

is_growable <- function(tree) return(length(tree$valid_grow_nodes) > 0)

is_changeable <- function(tree) return(length(tree$valid_change_nodes) > 0)

fixTreeParms <- function(tree){
  tree$node_idx   <- which(unlist(lapply(tree$nodes, length)) > 0)
  tree$num_nodes  <- length(tree$node_idx)
  tree$num_leaves <- sum(unlist(lapply(tree$nodes, \(n) n$leaf)) > 0)
  tree$num_nogs   <- sum(unlist(lapply(tree$nodes, \(n) n$nog)) > 0)
  tree$max_depth  <- max(unlist(lapply(tree$nodes, \(n) n$depth)))
  return(tree)
}
