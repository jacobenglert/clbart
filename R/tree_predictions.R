# Get predictions for x given a single tree
get_tree_predictions <- function(tree, x) {

  # Allocate storage for predictions
  predictions <- rep(NA_real_, nrow(x))

  # Define traversal function
  get_tree_predictions_recursive <- function(tree, indices) {

    # If the current node is a leaf
    if (tree$var[1] == -1) {

      # Store predictions
      predictions[indices] <<- tree$value[1]

      # Return 1 since this node has no further children
      return(1)
    }

    # Otherwise check which indices go left at the current node
    goesLeft <- x[indices, tree$var[1]] <= tree$value[1]

    # For indices that go left, move on to the left child of the current node
    # Need to skip over the root
    # headOfLeftBranch <- tree[-1,]
    # n_nodes.left <- get_tree_predictions_recursive(headOfLeftBranch, indices[goesLeft])

    leftNode <- get_children(tree$node[1])[1]
    leftBranchNodes <- sapply(tree$node, \(n) leftNode %in% c(n, get_ancestors(n)))
    leftBranch <- tree[leftBranchNodes,]
    get_tree_predictions_recursive(leftBranch, indices[goesLeft])

    # For indices that go right, move on to the right child of the current node
    # Need to skip over the root and all nodes on the left side of the tree
    # headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
    # n_nodes.right <- get_tree_predictions_recursive(headOfRightBranch, indices[!goesLeft])

    rightNode <- get_children(tree$node[1])[2]
    rightBranchNodes <- sapply(tree$node, \(n) rightNode %in% c(n, get_ancestors(n)))
    rightBranch <- tree[rightBranchNodes,]
    get_tree_predictions_recursive(rightBranch, indices[!goesLeft])

    # Return the total number of nodes including this one and its children
    # return(1 + n_nodes.left + n_nodes.right)
  }

  # Get predictions
  get_tree_predictions_recursive(tree, seq_len(nrow(x)))

  return(predictions)
}

# get_tree_predictions2 <- function (tree, x) {
#   leaves <- tree$node[tree$var == -1]
#   predictions <- rep(NA_real_, nrow(x))
#   for (l in leaves) predictions[get_mapped_observations(tree, l, w)] <- tree$value[tree$node == l]
#   return (predictions)
# }


# Get predictions for x given multiple trees
get_forest_predictions <- function(forest, x, add = TRUE) {

  # Get predictions for each forest in a tree
  predictions <- by(data = forest,
                    INDICES = forest[,c('tree')],
                    FUN = get_tree_predictions,
                    x = x)

  # Get the sum-of-trees predictions (default)
  if (add) predictions <- Reduce('+', predictions)

  return(predictions)
}

#' Get predictions for a set of predictors given posterior samples of forests
#'
#' @param forest_posterior List of forest data frames returned by the \code{clbart} function.
#' @param x Matrix of covariates to make predictions for.
#'
#' @return Returns a matrix of predictions with \code{nrow(x)} rows and \code{opts$num_save} columns, where opts$num_save is the number of posterior samples taken during model-fitting.
#' @export
#'
#' @examples
#' # fit is an object returned by clbart()
#' # w is a matrix of covariates used in the BART component
#' # pred_post <- get_forest_posterior_predictions(fit$forest, w)
get_forest_posterior_predictions <- function (forest_posterior, x) {

  # Get sum-of-trees predictions for each forest in the posterior
  predictions <- by(data = forest_posterior,
                    INDICES = forest_posterior[,c('sample')],
                    FUN = get_forest_predictions,
                    x = x)

  # Concatenate predictions by columns
  predictions <- do.call(cbind, predictions)
  colnames(predictions) <- NULL

  return(predictions)
}
