# This package provides functions that quantify infectious disease outbreaks using a branching process,
# a stochastic process in which each individual in generation n produces a random number of individuals
# in generation n+1, continuing for some number of generations or until there are no individuals
# remaining (stochastic extinction).

# The random number of next-generation individuals produced by each individual is drawn from the
# offspring distribution, a discrete probability distribution with non-negative range. To model
# infectious disease outbreaks, it is common to use a negative binomial offspring distribution,
# parameterized by the mean R and dispersion parameter k. This parameterization is equivalent to
# using mu = R and size = k in R's "NegBinomial", e.g. dnbinom(x, mu=R, size=k) would give the density,
# i.e. the probability of exactly x transmissions from one individual



#' Probability of y total transmission directly from x independent infected individuals
#'
#' @param x number of infected individuals in generation n
#' @param y number of total transmissions in generation n+1
#' @param R mean of negative binomial offspring distribution
#' @param k dispersion of negative binomial offspring distribution
#'
#' @examples
#' # With 5 individuals in this generation, what is the probability of
#' # 0 to 15 transmissions in the next generation?
#' p(5, 0:15, R=0.2, k=0.1)
#' @export
p <- function(x,y,R,k) exp(lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k))
