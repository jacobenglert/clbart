# Compute the CLR data log-likelihood for a given node
comp_loglik <- function(y, sc, stratum, max_win, na_locs, sum = TRUE){

  denoms <- NA * na_locs
  denoms[!na_locs] <- sc
  logLik <- sc[y == 1] - matrixStats::colLogSumExps(denoms, na.rm = TRUE)

  if(sum) return(sum(as.numeric(logLik)))
  else return(as.numeric(logLik))
}

# Compute the gradient of the CLR likelihood for a given node
comp_grd <- function(y, z, sc, stratum, max_win, na_locs){

  nums <- denoms <- NA * na_locs
  nums[!na_locs] <- z * exp(sc)
  denoms[!na_locs] <- exp(sc)

  num <- colSums(nums, na.rm = TRUE)
  denom <- colSums(denoms, na.rm = TRUE)

  grd <- z[y == 1] - num / denom

  return(sum(grd))
}

# Compute Fisher's information of the CLR likelihood for a given node
comp_fisher <- function(y, z, sc, stratum, max_win, na_locs){

  a <- b <- c <- NA * na_locs

  a[!na_locs] <- exp(sc)
  b[!na_locs] <- z * exp(sc)
  c[!na_locs] <- z^2 * exp(sc)

  a <- colSums(a, na.rm = TRUE)
  b <- colSums(b, na.rm = TRUE)
  c <- colSums(c, na.rm = TRUE)

  I <- ((a * c) - (b^2)) / (a^2)

  return(sum(I))
}

# comp_mv (compute the mean and variance for the proposal distribution)
comp_mv <- function(m_start = 0, y, z, stratum, sc1 = 0, lambda = 0, sigma2_mu, max_win, na_locs){

  m <- m_start

  # Use Fisher scoring to update m and v
  grd <- comp_grd(y, z, sc = sc1 + z * (m + lambda), stratum, max_win, na_locs) - (m / sigma2_mu)
  I <- comp_fisher(y, z, sc = sc1 + z * (m + lambda), stratum, max_win, na_locs) + (1 / sigma2_mu)

  i <- 0
  while(abs(grd) > sqrt(I)/10){
    m <- m + (grd / I)
    grd <- comp_grd(y, z, sc1 + z * (m + lambda), stratum, max_win, na_locs) - (m / sigma2_mu)
    I <- comp_fisher(y, z, sc1 + z * (m + lambda), stratum, max_win, na_locs)  + (1 / sigma2_mu)
    if(i == 10 | is.nan(I) | is.nan(grd)){
      m <- stats::coef(survival::clogit(y ~ offset(sc1) + offset(z * lambda) + z + survival::strata(stratum)))
      grd <- comp_grd(y, z, sc1 + z * (m + lambda), stratum, max_win, na_locs) - (m / sigma2_mu)
      I <- comp_fisher(y, z, sc1 + z * (m + lambda), stratum, max_win, na_locs)  + (1 / sigma2_mu)
    }
    i <- i + 1
  }

  # Return the mean and standard deviation of the proposal distribution
  return(list(m = m, v = I^(-1/2)))
}
