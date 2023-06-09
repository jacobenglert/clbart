% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict_forest}
\alias{predict_forest}
\title{Obtain predictions of heterogeneous log-odds ratios from a posterior forest
within a \code{clbart} model fit.}
\usage{
predict_forest(forest = NULL, new_data = NULL)
}
\arguments{
\item{forest}{list of trees (part of a \code{clbart} model posterior)}

\item{new_data}{data frame for which to obtain predictions (should match the
format of the moderating covarite data frame \code{w} in the initial \code{clbart()}
call)}
}
\value{
a numeric vector of predicted log-odds ratios
}
\description{
Obtain predictions of heterogeneous log-odds ratios from a posterior forest
within a \code{clbart} model fit.
}
\examples{
'INSERT EXAMPLE HERE'
}
