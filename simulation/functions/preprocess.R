######################################################
# MFLFM
######################################################
# clustering - data preprocessing
######################################################
# 2019.09.12
######################################################
# by DY
######################################################

# given containMULT & containCOV
if( (containCOV == 0) & (containMULT == 0) ){
  
  CASE <- 1
  
}else if( (containCOV == 1) & (containMULT == 0) ){
  
  CASE <- 2
  
}else if( (containCOV == 0) & (containMULT == 1) ){
  
  CASE <- 3
  
}else{
  
  CASE <- 4
}




if(containCAT == 0){
  
  q <- 3
  
  Zdat <- matrix(NA, q, n)
  
  Zdat[1,] <- c( rep(1, round(n/2) ), rep(2, n-round(n/2) ) )
  Zdat[2,] <- c( rep(1,round(n/3)), rep(2,round(n/3)), rep(3,n-2*round(n/3)) ) 
  Zdat[3,] <- c( rep(1,round(n/4)), rep(2,round(n/4)), rep(3,round(n/4)), rep(4,n-3*round(n/4)) ) 
}

qq <- nrow(Zdat)
qdims <- sapply(1:qq, function(x){ max(Zdat[x,]) })
qstar <- sum( qdims )
qind <- c(0,cumsum(qdims))
qind1 <- c(0,cumsum(qdims - 1))

Zdatstar  <- matrix(NA, qstar, n)
for(ww in 1:qq){
  
  tempZ <- matrix(NA, qdims[ww], n)
  for(ii in 1:n){
    
    tempvec <- rep(0, qdims[ww])
    tempvec[Zdat[ww,ii]] <- 1
    tempZ[,ii] <- tempvec
  }
  
  wwind <- (qind[ww]+1):qind[ww+1]
  Zdatstar[wwind,] <- tempZ
}





if( CASE == 3 ){
  
  # hyperparameter
  a_alpha <- b_alpha <- a_beta <- b_beta <- a_u <- b_u <- a_w <- b_w <- a_sigma <- b_sigma <- a_psi <- b_psi <- 1
  a_nu <- b_nu <- 1
  
  rhoA <- cc+2
  WA = diag(cc)
  
  # parameter
  r <- 8
  M <- 20
  
  tempK <- sum(Kcs)
  
  if(containCAT == 1){
    
    ZZ0 <- matrix( rbinom( (tempK+(qstar-qq))*r, 1, 1/2 ), (tempK+(qstar-qq)), r )
    VV0 <- matrix( rnorm( (tempK+(qstar-qq))*r ), (tempK+(qstar-qq)), r )
    Eta0 <- matrix( rnorm( n*r ), r, n )
    
    THETA0 <- (ZZ0 * VV0) %*% Eta0
    THETAcs0 <- list()
    for(ii in 1:cc){
      
      tempind <- (cumsum(c(0,Kcs))[ii]+1):(cumsum(c(0,Kcs))[ii+1])
      THETAcs0[[ii]] <- THETA0[tempind,]
    }
    Sigmas0 <- 1/rgamma(tempK+(qstar-qq), a_sigma, b_sigma)
    
  }else{
    
    ZZ0 <- matrix( rbinom( tempK*r, 1, 1/2 ), tempK, r )
    VV0 <- matrix( rnorm( tempK*r ), tempK, r )
    Eta0 <- matrix( rnorm( n*r ), r, n )
    
    THETA0 <- (ZZ0 * VV0) %*% Eta0
    THETAcs0 <- list()
    for(ii in 1:cc){
      
      tempind <- (cumsum(c(0,Kcs))[ii]+1):(cumsum(c(0,Kcs))[ii+1])
      THETAcs0[[ii]] <- THETA0[tempind,]
    }
    Sigmas0 <- 1/rgamma(tempK, a_sigma, b_sigma)
  }
  
  psics0 <- rep(1,cc)

  nu0 <- 1
  VVV0 <- rbeta(M-1, 1, nu0)
  temp <- VVV0 * c(1,cumprod(1-VVV0)[-(M-1)])
  PP0 <- c(temp, prod(1-VVV0) ); PP0 <- PP0/sum(PP0)
  LL0 <- sample(1:M, n, replace=TRUE, prob=PP0)
  
  wms0 <- 1/rgamma(r, a_w, b_w)
  mm0 <- rnorm(r)
  
  alpha0 <- 1
  beta0 <- 1
  ums0 <- 1/rgamma(r, a_u, b_u)
  
  Z_mu0 <- list()
  for(i in 1:M){ Z_mu0[[i]] <- as.vector( rmvnorm(1, mm0, diag(wms0) ) ) }
  Mus0 <- list(); for(i in 1:n){ Mus0[[i]] <- Z_mu0[[ LL0[i] ]] }
  
  llamb <- 5; ppi <- 1/2
  
  # main
  timer0 = proc.time()[3];
  
  A = diag(cc)
  
  containBETA <- 1
  containCOV <- 0 
  containMULT <- 1
  
  B <- Bcs
  THETA0 <- THETAcs0
  psi0 <- psics0
  
  THETA <- THETA0
  psi <- psi0
  ZZ <- ZZ0
  VV <- VV0
  Sigmas <- Sigmas0
  Eta <- Eta0
  LL <- LL0
  VVV <- VVV0
  PP <- PP0
  Z_mu <- Z_mu0
  Mus <- Mus0
  nu <- nu0
  wms <- wms0
  mm <- mm0
  alpha <- alpha0
  beta <- beta0
  ums <- ums0
  
}else if(CASE == 4){
  
  # hyperparameter
  a_alpha <- b_alpha <- a_beta <- b_beta <- a_u <- b_u <- a_w <- b_w <- a_sigma <- b_sigma <- a_psi <- b_psi <- 1
  a_nu <- b_nu <- 1

  rhoA <- cc+2
  WA = diag(cc)
  
  # parameter
  r <- 8
  M <- 20
  
  tempK <- sum(Kcs)
  
  if(containCAT == 1){
    
    ZZ0 <- matrix( rbinom( (tempK+p+(qstar-qq))*r, 1, 1/2 ), (tempK+p+(qstar-qq)), r )
    VV0 <- matrix( rnorm( (tempK+p+(qstar-qq))*r ), (tempK+p+(qstar-qq)), r )
    Eta0 <- matrix( rnorm( n*r ), r, n )
    
    THETA0 <- (ZZ0 * VV0) %*% Eta0
    THETAcs0 <- list()
    for(ii in 1:cc){
      
      tempind <- (cumsum(c(0,Kcs))[ii]+1):(cumsum(c(0,Kcs))[ii+1])
      THETAcs0[[ii]] <- THETA0[tempind,]
    }
    Sigmas0 <- 1/rgamma((tempK+p+(qstar-qq)), a_sigma, b_sigma)
    
  }else{
    
    ZZ0 <- matrix( rbinom( (tempK+p)*r, 1, 1/2 ), (tempK+p), r )
    VV0 <- matrix( rnorm( (tempK+p)*r ), (tempK+p), r )
    Eta0 <- matrix( rnorm( n*r ), r, n )
    
    THETA0 <- (ZZ0 * VV0) %*% Eta0
    THETAcs0 <- list()
    for(ii in 1:cc){
      
      tempind <- (cumsum(c(0,Kcs))[ii]+1):(cumsum(c(0,Kcs))[ii+1])
      THETAcs0[[ii]] <- THETA0[tempind,]
    }
    Sigmas0 <- 1/rgamma((tempK+p), a_sigma, b_sigma)
  }

  psics0 <- rep(1,cc)
  
  nu0 <- 1
  VVV0 <- rbeta(M-1, 1, nu0)
  temp <- VVV0 * c(1,cumprod(1-VVV0)[-(M-1)])
  PP0 <- c(temp, prod(1-VVV0) ); PP0 <- PP0/sum(PP0)
  LL0 <- sample(1:M, n, replace=TRUE, prob=PP0)
  
  wms0 <- 1/rgamma(r, a_w, b_w)
  mm0 <- rnorm(r)
  
  alpha0 <- 1
  beta0 <- 1
  ums0 <- 1/rgamma(r, a_u, b_u)
  
  Z_mu0 <- list()
  for(i in 1:M){ Z_mu0[[i]] <- as.vector( rmvnorm(1, mm0, diag(wms0) ) ) }
  Mus0 <- list(); for(i in 1:n){ Mus0[[i]] <- Z_mu0[[ LL0[i] ]] }
  
  llamb <- 4; ppi <- 1/2
  
  # main
  timer0 = proc.time()[3]
  
  containBETA <- 1
  containCOV <- 1
  containMULT <- 1
  
  A = diag(cc)
  
  B <- Bcs
  THETA0 <- THETAcs0
  psi0 <- psics0
  
  THETA <- THETA0
  psi <- psi0
  ZZ <- ZZ0
  VV <- VV0
  Sigmas <- Sigmas0
  Eta <- Eta0
  LL <- LL0
  VVV <- VVV0
  PP <- PP0
  Z_mu <- Z_mu0
  Mus <- Mus0
  nu <- nu0
  wms <- wms0
  mm <- mm0
  alpha <- alpha0
  beta <- beta0
  ums <- ums0
}



##########################
# data preprocessing
##########################
if(containBETA == 0){ beta <- 1 }

tol <- 10^(-50)
eps <- 10^(-6)

r <- nrow(Eta)
p <- nrow(Xdat)

n <- ifelse( containMULT == 1, ncol(THETA[[1]]), ncol(THETA) )
cc <- ifelse( containMULT == 1, length(B), 1 )

if( containMULT == 1 ){
  
  tempTHETA <- THETA[[1]]
  for(i in 2:cc){ tempTHETA <- rbind(tempTHETA, THETA[[i]]) }
  
  KKs <- sapply( THETA, nrow )
  tempK <- sum(KKs) #KK * cc
  kkind <- c(0,cumsum(KKs))
  
}else{
  
  tempTHETA <- THETA
  
  KKs <- nrow(THETA)
  tempK <- sum(KKs)
  kkind <- c(0,cumsum(KKs))
}

if( containCOV == 1 ){
  
  ZETA <- rbind( tempTHETA, Xdat )
  tempKp <- tempK + p
  
}else{
  
  ZETA <- tempTHETA
  tempKp <- tempK
}




if( containMULT == 1 ){
  
  Zmu <- matrix(0, r, M)
  for(i in 1:M){ Zmu[,i] <- Z_mu[[i]] }
  
  NOBS <- matrix( 0, length(nobs[[1]]), cc )
  YDAT <- matrix( 0, length(Ydat[[1]]), cc )
  BBB <- NULL
  for(ii in 1:cc){ BBB <- cbind(BBB, B[[ii]]) }
  for(cci in 1:cc){
    
    NOBS[,cci] <- nobs[[cci]]
    YDAT[,cci] <- Ydat[[cci]]
  }
  
  nobs <- NOBS
  Ydat <- YDAT
  B <- BBB
  
}else{
  
  Zmu <- matrix(0, r, M)
  for(i in 1:M){ Zmu[,i] <- Z_mu[[i]] }
  
  nobs <- matrix(nobs, ncol = 1)
  Ydat <- matrix(Ydat, ncol = 1)
}

Hp <- 0
for(i in 1:tempKp){ Hp <- Hp + 1/i }


if( containMULT == 1 ){
  
  psi0 <- 1
  
}else{
  
  psi0 <- psi
  psi <- c(1,1)
  A = diag(cc)
}




gis0 <- matrix(NA, qstar-qq, n)
for(ww in 1:qq){
  
  tempq <- qdims[ww]
  
  wwind <- (qind[ww]+1):qind[ww+1]
  tempz <- Zdatstar[wwind,]
  
  tempg <- matrix(NA, qdims[ww]-1, n)
  for(ii in 1:n){
    
    if(tempz[tempq,ii]==1){
      
      tempg[,ii] <- -1
      
    }else{
      
      tempg[,ii] <- 1
      tempg[which(tempz[,ii]==1),ii] <- 2
    }
  }
  
  wwind1 <- (qind1[ww]+1):qind1[ww+1]
  gis0[wwind1,] <- tempg
}
gis <- gis0


if( containCAT == 1 ){
  
  ZETA <- rbind( ZETA, gis )
  tempKpq <- tempKp + (qstar-qq)
  
  # ZZ0 <- ZZ <- rbind(ZZ, matrix( rbinom( (qstar-qq)*r, 1, 1/2 ), (qstar-qq), r ))
  # VV0 <- VV <- rbind(VV, matrix( rnorm( (qstar-qq)*r ), (qstar-qq), r ))
  # Sigmas0 <- Sigmas <- c( Sigmas, 1/rgamma((qstar-qq), a_sigma, b_sigma) )
  
}else{
  
  tempKpq <- tempKp
}








