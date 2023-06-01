resample <- function(x, ...){
  return(x[sample.int(length(x), ...)])
}

get_children <- function(id){
  return(c(id*2, id*2 + 1))
}

get_parent <- function(id){
  return(ifelse(id == 1, numeric(0), floor(id / 2)))
}

get_sibling <- function(id){
  c.ids <- get_children(get_parent(id))
  return(c.ids[which(c.ids != id)])
}

p_split <- function(d, alpha = 0.95, beta = 2){
  return(alpha * ((1 + d)^(-beta)))
}
