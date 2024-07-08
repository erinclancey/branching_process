#' Infectious disease outbreak quantification using branching processes
#'
#' This package provides functions that quantify infectious disease outbreaks
#' using a branching process, a stochastic process in which each individual in
#' generation n produces a random number of individuals in generation n+1,
#' continuing for some number of generations or until there are no individuals
#' remaining (stochastic extinction).
#'
#' The random number of next-generation individuals produced by each individual
#' is drawn from the offspring distribution, a discrete probability distribution
#' with non-negative range. To model infectious disease outbreaks, it is common
#' to use a negative binomial offspring distribution, parameterized by the mean
#' `R` and dispersion parameter `k`. This parameterization is equivalent to
#' using mu = R and size = k in R's "NegBinomial", e.g. dnbinom(x, mu=R, size=k)
#' would give the density,  i.e. the probability of exactly x transmissions from
#' one individual.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats uniroot
## usethis namespace: end
NULL
