get_node_map <- function (z, y, s, offset, ind, hypers, start_value = 0) {
  
  m <- start_value
  z <- z[ind, , drop = FALSE]
  y <- y[ind]
  s <- match(s[ind], unique(s[ind]))
  offset <- offset[ind]
  
  # Use Fisher scoring to update m and v
  grd <- clr_grd_cpp(m, z, y, s, offset) - ((m - hypers$mu_mu) / (hypers$sigma_mu^2))
  I <- clr_fish_cpp(m, z, y, s, offset) + (1 / hypers$sigma_mu^2)
  
  i <- 0
  while (TRUE) {
    
    m <- m + (grd / I)
    grd <- clr_grd_cpp(m, z, y, s, offset) - ((m - hypers$mu_mu) / (hypers$sigma_mu^2))
    I <- clr_fish_cpp(m, z, y, s, offset) + (1 / hypers$sigma_mu^2)
    
    # If algorithm is struggling to converge, try an alternative optimizer.
    # This may be needed on occasion if start_value is poor.
    if (i == 10 | is.nan(I) | is.infinite(I) | is.nan(grd) | is.infinite(grd)) {
      
      log_posterior <- function(value) {
        clr_loglik_cpp(value, z, y, s, offset) - (((value - hypers$mu_mu)^2) / (2 * hypers$sigma_mu^2))
      }
      
      m <- stats::optimize(log_posterior, lower = -10, upper = 10, maximum = TRUE)$maximum
      I <- clr_fish_cpp(m, z, y, s, offset) + (1 / hypers$sigma_mu^2)
      
      break
      
    }
    i <- i + 1
    
    if (abs(grd) > sqrt(I) / 10) break
    
  }
  
  # Return the mean and standard deviation of the proposal distribution
  return(list(m = m, v = I^(-1/2)))
  
}
