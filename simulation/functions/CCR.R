library(gtools)

CCR <- function(DPresult, original){
  
  
  n <- length(DPresult)
  
  if( length( unique(DPresult) ) >= length( unique(original) ) ){
    
    list1 <- unique(DPresult)
    ind1 <- DPresult
    
    list2 <- unique(original)
    ind2 <- original
    
  }else{
    
    list1 <- unique(original)
    ind1 <- original
    
    list2 <- unique(DPresult)
    ind2 <- DPresult
  }
  
  
  tempperm <- permutations(n=length(list1), r=length(list2), v=list1, repeats.allowed=F)
  temprate <- rep(NA, nrow(tempperm))
  for(j in 1:nrow(tempperm)){
    
    temporigin <- ind2
    for( t in 1:length(list2) ){ temporigin[ ind2 == list2[t] ] <- tempperm[j,t] }
    
    temprate[j] <- sum(temporigin == ind1)/n * 100
  }
  
  return(max(temprate))
}

CCRcomplete <- function(DPresult, original){
  
  n <- nrow(DPresult)
  rate <- rep(NA,n)
  for(i in 1:n){ rate[i] <- CCR( DPresult[i,], original ) }
  
  return(rate)
}