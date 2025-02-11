
<!-- README.md is generated from README.Rmd. Please edit that file -->

# clbart

The goal of `clbart` is to estimate heterogeneous effects of a primary
exposure in a scenario where only cases are observed. Specifically, this
method is appropriate for data which have been constructed using the
design, such that only one case exists within each strata. The effect of
the primary exposure is modeled as a function of time-invariant
covariates using Bayesian additive regression trees.

## Installation

You can install the development version of `clbart` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jacobenglert/clbart")
```

## Example

For a detailed example of how to use the package and fit `clbart`
models, view the package vignette
[here](https://github.com/jacobenglert/clbart/blob/main/vignettes/clbart.pdf).
