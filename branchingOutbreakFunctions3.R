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

p <- function(x,y,R,k) exp(lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k))

# Example:
#
# p(5, 0:15, R=0.2, k=0.1)
#
# With 5 individuals in generation n and negative binomial offspring distribution with mean R=0.2
# and dispersion k=0.1, gives the probability that generation n+1 has a total of 0 to 15 individuals 

# qAny is the probability that n initial cases lead to an extinguished outbreak of total size j
# after any number of transmission generations (j includes the n initial cases)

qAny <- function(n,j,R,k) p(j,j-n,R,k)*n/j

# Example:
# 
# q(5, 5:20, R=0.2, k=0.1)
#
# With 5 initial individuals and negative binomial offspring distribution with mean R=0.2
# and dispersion k=0.1, gives the probability of outbreak extinction with a total number
# final outbreak size of exactly 5 to 20 inidividuals (including the initial 5). 


# pgl[g] is the probability that one initial case leads to an outbreak lasting less than g generations of transmission
# getPgl calculates pgl[g] for all g from 1 to gMax

getPgl <- function(gMax,R,k){
  pgl <- rep(0,gMax)
  pgl[1] <- (1+R/k)^(-k)
  if(gMax > 1) for(g in 2:gMax) pgl[g] <- (1+R/k*(1-pgl[g-1]))^(-k)	
  pgl
}

# pExtinct is the outbreak extinction probability from one initial case
# pgl[g] will converge to pExtinct for large enough g

pExtinct <- function(R,k) ifelse(R > 1, uniroot(function(q) q - (1+R/k*(1-q))^(-k), c(0,0.999999), tol=1e-7)$root, 1)

# qg is the probability that n initial cases lead to an outbreak that ends after
# exactly g generations of transmission AND has exactly j (> n+g-1) total cases 

qg <- function(g,n,j,R,k){
  
  if(g==1){
    out <- p(n,j-n,R,k)*p(j-n,0,R,k)
  }else if(g==2){
    out <- sum(p(n,1:(j-n-1),R,k) * p(1:(j-n-1),(j-n-1):1,R,k) * p((j-n-1):1,0,R,k))
  }else{
    
    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))
    
    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1
    
    pProd <- p(n,x1,R,k)
    
    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd * p(x[,i-1],x[,i],R,k)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * p(x[,g-1],xLast,R,k) * p(xLast,0,R,k)
    out <- sum(pProd)
  }
  out
}

# pg is the probability that n initial cases lead to an outbreak that lasts at least g generations
# of transmission AND has exactly j (> n+g-1) total cases after generation g

pg <- function(g,n,j,R,k){
  
  if(g==1){
    out <- p(n,j-n,R,k)
  }else if(g==2){
    out <- sum(p(n,1:(j-n-1),R,k) * p(1:(j-n-1),(j-n-1):1,R,k))
  }else{
    
    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))
    
    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1
    
    pProd <- p(n,x1,R,k)
    
    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd * p(x[,i-1],x[,i],R,k)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * p(x[,g-1],xLast,R,k)
    out <- sum(pProd)
  }
  out
}


# The following are new versions of some of the above formulas, for a model in which (R,k) = (R0,k0) for the
# initial case(s) ONLY, while any cases in subsequent generations transmit according to (R,k) = (Rc,kc)

qAnySwitch1 <- function(n,j,R0,k0,Rc,kc) ifelse(j==n, p(n,0,R0,k0), sum(p(n,1:(j-n),R0,k0) * qAny(1:(j-n),j-n,Rc,kc)))

getPglSwitch1 <- function(gMax,R0,k0,Rc,kc){
  pgl <- rep(0,gMax)
  pgl[1] <- (1+R0/k0)^(-k0)
  if(gMax > 1) for(g in 2:gMax) pgl[g] <- (1+Rc/kc*(1-pgl[g-1]))^(-kc)	
  pgl
}

qgSwitch1 <- function(g,n,j,R0,k0,Rc,kc){
  
  if(g==1){
    out <- p(n,j-n,R0,k0)*p(j-n,0,Rc,kc)
  }else if(g==2){
    out <- sum(p(n,1:(j-n-1),R0,k0) * p(1:(j-n-1),(j-n-1):1,Rc,kc) * p((j-n-1):1,0,Rc,kc))
  }else{
    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))
    
    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1
    
    pProd <- p(n,x1,R0,k0)
    
    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd * p(x[,i-1],x[,i],Rc,kc)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * p(x[,g-1],xLast,Rc,kc) * p(xLast,0,Rc,kc)
    out <- sum(pProd)
  }
  out
}

pgSwitch1 <- function(g,n,j,R0,k0,Rc,kc){
  
  if(g==1){
    out <- p(n,j-n,R0,k0)
  }else if(g==2){
    out <- sum(p(n,1:(j-n-1),R0,k0) * p(1:(j-n-1),(j-n-1):1,Rc,kc))
  }else{
    
    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))
    
    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1
    
    pProd <- p(n,x1,R0,k0)
    
    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd * p(x[,i-1],x[,i],Rc,kc)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * p(x[,g-1],xLast,Rc,kc)
    out <- sum(pProd)
  }
  out
}

