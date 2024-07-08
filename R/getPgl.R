# pgl[g] is the probability that one initial case leads to an outbreak lasting less than g generations of transmission
# getPgl calculates pgl[g] for all g from 1 to gMax

getPgl <- function(gMax,R,k){
  pgl <- rep(0,gMax)
  pgl[1] <- (1+R/k)^(-k)
  if(gMax > 1) for(g in 2:gMax) pgl[g] <- (1+R/k*(1-pgl[g-1]))^(-k)	
  pgl
}
