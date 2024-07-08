#' Probability of outbreak extinction
#'
#' qAny is the probability that n initial cases lead to an extinguished outbreak
#' of total size j after any number of transmission generations (j includes the
#' n initial cases)
#'
#' @param n number of initial cases in generation 0
#' @param j total outbreak size (>= n).
#' @param R mean of negative binomial offspring distribution
#' @param k dispersion of negative binomial offspring distribution
#' @examples
#' # With 5 initial individuals and negative binomial offspring distribution with mean R=0.2
#' # and dispersion k=0.1, gives the probability of outbreak extinction with a total number
#' # final outbreak size of exactly 5 to 20 inidividuals (including the initial 5).
#' qAny(5, 5:20, R=0.2, k=0.1)
#'
#' @export
qAny <- function(n,j,R,k) p(j,j-n,R,k)*n/j
