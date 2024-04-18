#' Predict Method for CL-BART Model Fits
#'
#' @param object a fitted object of class \code{clbart}.
#' @param newdata optionally, a list of data frames to be used for prediction.
#' When \code{type = "con"} is requested, this list must include a data frame of
#' covariates \code{x}. When \code{type = "bart"} is requested, this list must
#' include a data frame of exposure moderating covariates \code{w}.
#' @param type the type of prediction requested. The default (\code{"bart"})
#' produces predictions from the BART component of the CL-BART model fit on
#' the log odds-ratio scale, and \code{type = "con"} returns estimates of the
#' parametric component of the CL-BART model fit on the log odds scale.
#' @param posterior logical indicating whether or not the entire posterior
#' distribution for the requested predictions be returned. if \code{TRUE}
#' (default), a matrix is returned where each column represents the predictions
#' for a single posterior sample. If \code{FALSE}, a vector of posterior means
#' is returned.
#'
#' @return a vector (or matrix if \code{posterior = TRUE}) of predictions from
#' a fitted \code{clbart} model object.
#' @export
#'
predict.clbart <- function (object,
                            newdata,
                            type = c('bart','con'),
                            posterior = FALSE) {

  if (type == 'bart') {
    if (is.null(newdata$w)) stop("w data frame must be specified in newdata.")
    if (is.null(fit$forest)) stop("object is missing forest component.")
    if (ncol(fit$split_counts) == ncol(newdata$w)) {
      if (any(colnames(fit$split_counts) != colnames(newdata$w))) {
        stop("Moderating covariates in newdata$w do not match those used in fitted clbart model.")
      }
    } else {
      stop("Number of columns in newdata$w does not match number of moderators in fitted clbart model.")
    }
    predictions <- get_forest_posterior_predictions(object$forest, newdata$w)
  } else if (type == 'con') {
    if (is.null(newdata$x)) stop("x data frame must be specified in newdata.")
    if (ncol(fit$beta) == ncol(newdata$x)) {
      if (any(colnames(fit$beta) != colnames(newdata$x))) {
        stop("Confounding covariates in newdata$x do not match those in object$beta.")
      }
    } else {
      stop("Number of columns in newdata$x does not match number of columns in object$beta.")
    }
    predictions <- as.matrix(newdata$x) %*% t(object$beta)
  } else {
    stop("Unsupported type specified.")
  }

  if (!posterior) predictions <- rowMeans(predictions)

  return(predictions)
}

