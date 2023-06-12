
predict_tree <- function(tree = NULL, new_data = NULL){
  leaves <- which(lapply(tree$nodes, \(n) as.numeric(n$leaf)) == 1)
  if(is.null(new_data)){
    stopifnot('Must supply new data for prediction' = !is.null(tree$nodes[[1]]$data_idx))
    return(Reduce('+', lapply(tree$nodes[leaves], \(l) l$mu * l$data_idx)))
  }
  else{
    pred <- numeric(nrow(new_data))

    if(length(leaves) == 1){
      pred[1:length(pred)] <- tree$nodes[[1]]$mu
    }
    else{
      design_mat <- lapply(tree$nodes[leaves], \(l){
        eval(parse(text = paste0('new_data$', l$split_rule, collapse = ' & ')))}) |>
        do.call(what = cbind)
      mu_vec <- lapply(tree$nodes[leaves], '[[', 'mu') |> do.call(what = rbind)

      pred <- design_mat %*% mu_vec |> as.numeric()
    }
    return(pred)
  }
}

#' Obtain predictions of heterogeneous log-odds ratios from a posterior forest
#' within a `clbart` model fit.
#'
#' @param forest list of trees (part of a `clbart` model posterior)
#' @param new_data data frame for which to obtain predictions (should match the
#' format of the moderating covarite data frame `w` in the initial `clbart()`
#' call)
#'
#' @return a numeric vector of predicted log-odds ratios
#' @export
#'
#' @examples
#' 'INSERT EXAMPLE HERE'
predict_forest <- function(forest = NULL, new_data = NULL){
  return(Reduce('+', lapply(forest, \(t) predict_tree(t, new_data))))
}

simplify_forest <- function(forest){

  new_forest <- list()

  for(i in 1:length(forest)){
    new_forest[[i]] <- list(nodes = lapply(forest[[i]]$nodes, \(n) n[c('leaf','split_rule','mu')]))
  }

  return(new_forest)
}
