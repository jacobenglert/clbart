get_mapped_observations <- function (tree, node, x) {
  
  n <- nrow(x)

  if (node == 1) return(rep(TRUE, n))
    
  # Retrieve branch
  anc <- get_ancestors(node)
  branch <- tree[tree$node %in% anc,]
  wentLeft <- c(anc[-1], node) %% 2 == 0
  
  # Iterate through splitting rules
  indices <- 1:n
  for (i in 1:nrow(branch)) {
    if (wentLeft[i]) {
      indices <- indices[which(x[indices, branch$var[i]] <= branch$value[i])]
    } else {
      indices <- indices[which(x[indices, branch$var[i]] > branch$value[i])]
    }
  }
  
  is_mapped <- logical(n)
  is_mapped[indices] <- TRUE
  
  return(is_mapped)
}
  