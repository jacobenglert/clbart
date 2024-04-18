#' Summary Method for CL-BART Fits
#'
#' @description
#' Produces basic summary statistics obtained from a fitted \code{clbart} model.
#'
#'
#' @param object a fitted object of class \code{clbart}.
#' @param ... additional arguments affecting the summary produced.
#'   - `digits`: integer, specifying the number of decimal places (default is 2).
#'
#' @return A list of basic summaries for different components of a fitted
#' \code{clbart} model. The summary includes:
#' - posterior means and credible intervals for traditional regression
#' coefficients (on the log odds-ratio scale).
#' - posterior mean and credible interval for the overal average exposure effect
#' (on the log odds-ratio scale).
#' - variable importance measures including posterior mean split proportions,
#' posterior mean splitting probabilities, and posterior inclusion probabilities
#'  for each exposure effect moderator.
#' - additional information (e.g., WAIC).
#' @export
#' @examples
#' \dontrun{
#' # Assuming `model` is a fitted model of class `clbart`
#' summary(model, digits = 3)
#' }
summary.clbart <- function (object, ...) {

  # Handle additional arguments
  args <- list(...)
  digits <- if("digits" %in% names(args)) args$digits else 2

  if(!is.null(object$beta)) {
    beta_stats <- data.frame(
      post.mean = round(colMeans(object$beta), digits),
      post.2.5 = round(apply(object$beta, 2, \(x) stats::quantile(x, 0.025)), digits),
      post.97.5 = round(apply(object$beta, 2, \(x) stats::quantile(x, 0.975)), digits)
    )
    rownames(beta_stats) <- colnames(object$beta)
  } else {
    beta_stats <- NULL
  }

  bart_stats <- data.frame(
    post.mean <- round(mean(object$lambda_mean_overall), digits),
    post.2.5 <- round(stats::quantile(object$lambda_mean_overall, 0.025), digits),
    post.97.5 <- round(stats::quantile(object$lambda_mean_overall, 0.975), digits)
  )
  rownames(bart_stats) <- 'average bart prediction'
  colnames(bart_stats) <- c('post.mean','post.2.5','post.97.5')

  splits <- object$split_counts / rowSums(object$split_counts)
  splits[is.nan(splits)] <- 0
  var_imp <- data.frame(

    split.props <- round(colMeans(splits), digits),
    split.probs <- round(colMeans(object$split_probs), digits),
    var.pip <- round(colMeans(object$split_counts > 0), digits)
  )
  rownames(var_imp) <- colnames(object$split_counts)
  colnames(var_imp) <- c('split.props','split.probs','PiP')

  summary <- list(
    beta_stats = beta_stats,
    bart_stats = bart_stats,
    var_imp = var_imp,

    WAIC = round(object$WAIC),
    pWAIC = round(object$pWAIC)
  )

  class(summary) <- 'summary_clbart'

  return(summary)
}

print.summary_clbart <- function(x, ...) {
  cat("Summary of CL-BART Model:\n\n")

  cat("Confounder Coefficients:\n")
  print(x$beta_stats)

  cat("\nBART Predictions:\n")
  print(x$bart_stats)

  cat("\nBART Variable Importance:\n")
  print(x$var_imp)

  cat("\nOther Information:\n")
  print(paste0('WAIC: ', x$WAIC, '; pWAIC: ', x$pWAIC))

  invisible(x)
}
