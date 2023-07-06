#' Simulated Case-Crossover Data
#'
#' A simulated dataset representing a case-crossover analysis with 100 subjects
#' observed over 5 time points each.
#'
#' @format ## `cco`
#' A data frame with 500 rows and 21 columns:
#' \describe{
#'   \item{X1, ..., X9}{Time-varying covariates.}
#'   \item{W1, ..., W9}{Within subject time-invariant covariates}
#'   \item{Y}{Binary observed outcome. Only one event per subject due to the 1:5
#'   case-crossover matching scheme.}
#'   \item{Z}{Time-varying exposure variable.}
#'   \item{Strata}{Identifier for subject.}
#' }
"cco"
