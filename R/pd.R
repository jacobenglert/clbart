
#' Partial dependence for clbart
#'
#' Compute partial dependence statistics
#'
#' @param x Data frame of covariates used in the model-fitting process.
#' @param f_hat Function that returns a posterior distribution of predictions given covariates \code{x}.
#' @param vars Character vector of covariate names to vary in partial dependence calculations.
#' @param K Numeric vector of unique values to evaluate each covariate in \code{vars} at. By default this value will be decreased to the number of unique values if necessary.
#' @param f (Optional) function that returns the true prediction for a set of covariates. Useful when conducting simulation studies.
#'
#' @return A data frame of posterior mean estimates and 95% credible intervals.
#' @export
#'
#' @examples
#' # fit is an object returned by clbart()
#' # w is a matrix of covariates used in the BART component during model-fitting
#' # f_hat <- function (x) get_forest_posterior_predictions(fit$forest, x)
#' # pdw1 <- pd(w, f_hat, vars = 'W1', K = 10)

pd <- function (x, f_hat, vars, K, f = NULL) {

  x <- as.data.frame(x)
  p <- length(vars)

  if (K %% 1 != 0) stop("K must be a positive integer.")
  if (length(K) == 1) K <- rep(K, p)

  # Create grid of variables to explore
  var.grids <- list()
  for (j in seq_len(p)) {
    var <- vars[j]
    v <- x[, var, drop = TRUE]
    n_unique <- min(length(unique(v)), K[j])
    var.grids[[var]] <- seq(min(v), max(v), length.out = n_unique)
  }
  grid <- expand.grid(var.grids)

  # Calculate partial dependence for each grid cell
  for (row in seq_len(nrow(grid))) {
    pdata <- x
    for (var in vars) pdata[[var]] <- grid[[var]][row]
    preds <- colMeans(f_hat(pdata))
    grid$est[row] <- mean(preds)
    grid$lcl[row] <- stats::quantile(preds, 0.025)
    grid$ucl[row] <- stats::quantile(preds, 0.975)

    # Calculate the true partial dependence if provided
    if (!is.null(f)) grid$truth[row] <- mean(f(pdata))
  }

  return(grid)

}
