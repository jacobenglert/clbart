cl_loglik <- function(beta, y, x, strata, na_locs, offset = numeric(nrow(x))){
  a <- NA * na_locs
  a[!na_locs] <- exp(x %*% beta + offset)
  loglik <- sum(x[y == 1,] %*% beta + offset[y == 1] - log(colSums(a, na.rm = TRUE)))
  return(loglik)
}

cl_grd <- function(beta, y, x, strata, na_locs, offset = numeric(nrow(x))){
  exbeta <- exp(x %*% beta + offset)
  b <- NA * na_locs
  b[!na_locs] <- exbeta
  b1 <- colSums(b, na.rm = TRUE)
  grd <- numeric(length(beta))
  for(j in 1:length(beta)){
    a <- NA * na_locs
    a[!na_locs] <- (x[,j]) * exbeta
    grd[j] <- sum(x[y == 1, j] - colSums(a, na.rm = TRUE) / b1)
  }
  return(grd)
}

cl_fish <- function(beta, x, strata, na_locs, offset = numeric(nrow(x))){
  exbeta <- exp(x %*% beta + offset)
  a <- NA * na_locs
  a[!na_locs] <- exbeta
  a1 <- colSums(a, na.rm = TRUE)

  I <- matrix(0, ncol = length(beta), nrow = length(beta))
  for(j in 1:length(beta)){
    c <- NA * na_locs
    c[!na_locs] <- x[,j] * exbeta
    c1 <- colSums(c, na.rm = TRUE)

    for(j2 in 1:length(beta)){
      b <- d <- NA * na_locs
      b[!na_locs] <- x[,j] * x[,j2] * exbeta
      d[!na_locs] <- x[,j2] * exbeta

      I[j,j2] <- sum(((a1 *  colSums(b, na.rm = TRUE)) - (c1 *  colSums(d, na.rm = TRUE))) / a1^2)
    }
  }
  return(I)
}

clr <- function(y, x, strata, offset = numeric(nrow(x))){

  # Identify strata lengths for faster computation
  strata <- as.numeric(as.factor(strata))
  n <- length(unique(strata))
  windows <- table(strata)
  max_win <- max(windows)
  na_locs <- (1:(max_win * n)) %in% (max_win * unique(strata)[which(windows != max_win)]) |>
    matrix(nrow = max_win, ncol = n)

  # Initialize beta vector, log-likelihood, gradient, and Fisher information
  beta <- matrix(0, ncol = 1, nrow = ncol(x))
  loglik <- loglik_old <- cl_loglik(beta, y, x, strata, na_locs, offset)
  grd <- cl_grd(beta, y, x, strata, na_locs, offset)
  I <- cl_fish(beta, x, strata, na_locs, offset)

  # Iterate until convergence
  s <- 1 # initial step size
  while(abs(loglik - loglik_old) / abs(loglik_old) < 1e-7){

    beta_old <- beta
    loglik_old <- loglik

    beta <- beta + s * solve(I) %*% grd
    loglik <- cl_loglik(beta, y, x, strata, na_locs, offset)
    grd <- cl_grd(beta, y, x, strata, na_locs, offset)
    I <- cl_fish(beta, x, strata, na_locs, offset)

    s <- s / 2 # half the step size
  }

  # Export estimated coefficients and variance-covariance matrix
  colnames(I) <- rownames(I) <- colnames(x)
  beta <- as.numeric(beta)
  names(beta) <- colnames(x)
  return(list(beta = beta, vcov = solve(I)))
}
