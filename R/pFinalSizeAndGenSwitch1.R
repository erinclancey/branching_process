#' Joint probability of outbreak final size and number of transmission generations
#'
#' @param g number of generations
#' @param n number of initial cases
#' @param j final size
#' @param R0 mean of negative binomial offspring distribution from generation one
#' @param k0 dispersion of negative binomial offspring distribution from generation one
#' @param Rc mean of negative binomial offspring distribution from generation two on
#' @param kc dispersion of negative binomial offspring distribution from generation two on
#' @author Damon Toth
#' @export
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