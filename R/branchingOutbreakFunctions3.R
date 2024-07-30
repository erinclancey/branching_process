

# pExtinct is the outbreak extinction probability from one initial case
# pgl[g] will converge to pExtinct for large enough g

pExtinct <- function(R,k) ifelse(R > 1, uniroot(function(q) q - (1+R/k*(1-q))^(-k), c(0,0.999999), tol=1e-7)$root, 1)

# qg is the probability that n initial cases lead to an outbreak that ends after
# exactly g generations of transmission AND has exactly j (> n+g-1) total cases

pFinalSizeAndGen <- function(g,n,j,R,k){

  if(g==1){
    out <- pNextGen(n,j-n,R,k)*pNextGen(j-n,0,R,k)
  }else if(g==2){
    out <- sum(pNextGen(n,1:(j-n-1),R,k) * pNextGen(1:(j-n-1),(j-n-1):1,R,k) * pNextGen((j-n-1):1,0,R,k))
  }else{

    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))

    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1

    pProd <- pNextGen(n,x1,R,k)

    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd * pNextGen(x[,i-1],x[,i],R,k)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * pNextGen(x[,g-1],xLast,R,k) * pNextGen(xLast,0,R,k)
    out <- sum(pProd)
  }
  out
}

# pg is the probability that n initial cases lead to an outbreak that lasts at least g generations
# of transmission AND has exactly j (> n+g-1) total cases after generation g

pSizeAtGen <- function(g,n,j,R,k){

  if(g==1){
    out <- pNextGen(n,j-n,R,k)
  }else if(g==2){
    out <- sum(pNextGen(n,1:(j-n-1),R,k) * pNextGen(1:(j-n-1),(j-n-1):1,R,k))
  }else{

    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))

    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1

    pProd <- pNextGen(n,x1,R,k)

    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd * pNextGen(x[,i-1],x[,i],R,k)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd * pNextGen(x[,g-1],xLast,R,k)
    out <- sum(pProd)
  }
  out
}



# The following are new versions of some of the above formulas, for a model in which (R,k) = (R0,k0) for the
# initial case(s) ONLY, while any cases in subsequent generations transmit according to (R,k) = (Rc,kc)

pFinalSizeSwitch1 <- function(n,j,R0,k0,Rc,kc)
  ifelse(j==n, pNextGen(n,0,R0,k0), sum(pNextGen(n,1:(j-n),R0,k0) * pFinalSize(1:(j-n),j-n,Rc,kc)))

pGenSwitch1 <- function(gMax,R0,k0,Rc,kc){
  pgl <- rep(0,gMax)
  pgl[1] <- (1+R0/k0)^(-k0)
  if(gMax > 1) for(g in 2:gMax) pgl[g] <- (1+Rc/kc*(1-pgl[g-1]))^(-kc)
  pgl
}

pFinalSizeAndGenSwitch1 <- function(g,n,j,R0,k0,Rc,kc){

  if(g==1){
    out <- pNextGen(n,j-n,R0,k0)*pNextGen(j-n,0,Rc,kc)
  }else if(g==2){
    out <- sum(pNextGen(n,1:(j-n-1),R0,k0) * pNextGen(1:(j-n-1),(j-n-1):1,Rc,kc) * pNextGen((j-n-1):1,0,Rc,kc))
  }else{
    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))

    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1

    pProd <-pNextGen(n,x1,R0,k0)

    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd *pNextGen(x[,i-1],x[,i],Rc,kc)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd *pNextGen(x[,g-1],xLast,Rc,kc) *pNextGen(xLast,0,Rc,kc)
    out <- sum(pProd)
  }
  out
}

pSizeAtGenSwitch1 <- function(g,n,j,R0,k0,Rc,kc){

  if(g==1){
    out <-pNextGen(n,j-n,R0,k0)
  }else if(g==2){
    out <- sum(pNextGen(n,1:(j-n-1),R0,k0) *pNextGen(1:(j-n-1),(j-n-1):1,Rc,kc))
  }else{

    rs1 <- (j-n-g+1):1
    x1 <- rep(1:(j-n-g+1),choose(rs1+g-3,g-2))

    x <- matrix(0,length(x1),g-1)
    x[,1] <- x1

    pProd <-pNextGen(n,x1,R0,k0)

    rsA <- rs1
    for(i in 2:(g-1)){
      rsB <- sequence(rsA,rsA,-1)
      x[,i] <- rep(sequence(rsA),choose(rsB+g-2-i,g-1-i))
      pProd <- pProd *pNextGen(x[,i-1],x[,i],Rc,kc)
      rsA <- rsB
    }
    xLast <- j-n-rowSums(x)
    pProd <- pProd *pNextGen(x[,g-1],xLast,Rc,kc)
    out <- sum(pProd)
  }
  out
}

