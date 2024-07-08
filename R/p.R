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

# p is the probability of y total transmission directly from x independent infected individuals:



# Example:
#
# p(5, 0:15, R=0.2, k=0.1)
#
# With 5 individuals in generation n and negative binomial offspring distribution with mean R=0.2
# and dispersion k=0.1, gives the probability that generation n+1 has a total of 0 to 15 individuals 

p <- function(x,y,R,k) exp(lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k))