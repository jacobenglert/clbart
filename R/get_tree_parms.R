get_split_props <- function(forest, w){

  split_freqs <- table(factor(
    unlist(lapply(forest, \(t) lapply(t$nodes, '[[', 'split_var'))),
    levels = colnames(w)))
  split_props <- split_freqs / sum(split_freqs)

  return(split_props)

}

get_tree_props <- function(forest, w){

  tree_props <- lapply(forest, \(t) table(factor(
    unlist(lapply(t$nodes, '[[', 'split_var')), levels = colnames(w))) > 1) |>
    do.call(what = rbind) |>
    colMeans()

  return(tree_props)

}

get_avg_tree_depth <- function(forest){
  return(mean(unlist(lapply(forest, '[[', 'max_depth'))))
}

get_avg_num_nodes <- function(forest){
  return(mean(unlist(lapply(forest, '[[', 'num_nodes'))))
}

get_avg_num_leaves <- function(forest){
  return(mean(unlist(lapply(forest, '[[', 'num_leaves'))))
}
