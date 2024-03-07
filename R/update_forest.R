
update_forest <- function (forest, w, y, z, strata, offset, hypers) {
  
  x_beta <- offset
  lambda <- get_forest_predictions(forest, w)
  
  # old_forest <- forest
  # while (forest$var[1] != 1) {
  #   forest <- grow(old_forest, w, y, z, strata, x_beta, hypers)
  # }
  # old_forest <- forest
  # while (forest$var[2] != 2) {
  #   forest <- grow(old_forest, w, y, z, strata, x_beta, hypers)
  # }
  # old_forest <- forest
  # while (forest$var[3] != 3) {
  #   forest <- grow(old_forest, w, y, z, strata, x_beta, hypers)
  # }
  
  # Iterate through trees
  for (t in seq_len(hypers$num_tree)) {

    # Horseshoe prior adjustment
    # hypers$sigma_mu <- hypers$omega * hypers$tau[t]
    
    # Retrieve the current tree
    curr_tree <- forest[forest$tree == t,]
    
    # Subtract the prediction of the current tree
    lambda_r <- lambda - get_tree_predictions(curr_tree, w)
    
    # Propose a new tree
    move <- sample(1:3, size = 1, prob = p_moves(curr_tree))

    if (move == 1) {
      curr_tree <- grow(curr_tree, w, y, z, strata, x_beta + z * lambda_r, hypers)
    } else if (move == 2) {
      curr_tree <- prune(curr_tree, w, y, z, strata, x_beta + z * lambda_r, hypers)
    } else if (move == 3) {
      curr_tree <- change(curr_tree, w, y, z, strata, x_beta + z * lambda_r, hypers)
    }
    
    # Sample terminal node values (regardless of whether tree update was accepted)
    curr_tree <- update_mu(curr_tree, w, y, z, strata, x_beta + z * lambda_r, hypers)
    
    # Replace the previous tree with the updated tree
    forest <- rbind(forest[forest$tree < t,], curr_tree, forest[forest$tree > t,])
    
    # Add back prediction
    lambda <- lambda_r + get_tree_predictions(curr_tree, w)
  }
  
  return(forest)
  
}
