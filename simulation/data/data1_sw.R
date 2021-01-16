library(mvtnorm)
library(fda); library(KFAS); library(MCMCpack); 
#library(stochvol);
library(truncnorm) ; 
library(dlm)

getSplineInfo = function(tau, KNOTS=20, intKnots = NULL){
  
  allTaus = sort(unique(tau)) 	# all observation points
  a = min(allTaus)        	# lower endpoint
  b = max(allTaus)        	# upper endpoint
  
  numIntKnots = KNOTS 	# number of interior knots (= M in paper)
  if( is.null(intKnots) ){ intKnots = quantile(allTaus,seq(0,1,length= (numIntKnots+2))[-c(1,(numIntKnots+2))])}  	# interior knots 
  
  basis = create.bspline.basis(c(a,b),breaks=c(a,intKnots,b))						# spline basis
  blin = create.monomial.basis(c(a,b), nbasis=2) 								# linear basis
  
  Phi = bs(allTaus,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE) 		# basis matrix
  Omega = eval.penalty(basis, Lfdobj=2, rng=c(a,b)) 					 		# derivative matrix
  Jkl = eval.penalty(basis, Lfdobj=0, rng=c(a,b))						 		# integral matrix
  
  # Now, transform these matrices to align with the linear/nonlinear decomposition of d_k:
  # The first two entries of d_k are the linear components, the remaining are nonlinear (see Wand and Ormerod, 2008)
  eigOmega = eigen(Omega)
  indsZ = 1:(numIntKnots+2)
  UZ = eigOmega$vectors[, indsZ] 			# Linear part
  LZ = t(t(UZ)/sqrt(eigOmega$values[indsZ])) 	# Nonlinear part
  
  # Basis matrices
  PhiL = cbind(1, allTaus)				# Linear part
  PhiN = Phi%*%LZ						# Nonlinear part
  Phi = cbind(PhiL, PhiN)
  
  return(Phi)
}

# Preda
n1 <- round(n/5)
n2 <- round(n/5)
n3 <- round(n/5)
n4 <- round(n/5)
n5 <- n-n1-n2-n3-n4
te <- seq(0,1,length.out = ss) #seq(1,21,0.05)

# observation numbers
nobs1 <- nobs2 <- nobs3 <- rep(length(te),n) #sample(40:60, n, replace=TRUE)
nobs <- list(nobs1, nobs2, nobs3)

p <- 8
Xdat <- matrix(NA, p, n)

q <-3
Zdat <- matrix(NA, q, n)


# fval1 <- 2*log( te + 1/10 )
# fval2 <- 4*cos(5*te)+2
# fval3 <- exp(te)
# fval4 <- tan(te+1/4)
# fval5 <- 5*te^3 + 1
# fval6 <-3*sin( 3*te + 1/3 ) - 1


fval1 <- 2*log( te + 1/10 )+2
fval2 <- 5/2*cos(5*te)+2 - 3/2
fval3 <- exp(te)
fval4 <- tan(te+1/3)
fval5 <- 5*sin( 3*te + 1/2 ) - 2
fval6 <- te^3 - 3*te^2



fval7 <- te^2 - 5*te^4
fval8 <- exp(2*te^2)
fval9 <- 5*cos( 5*te^3+4.5 )


# fval1 <- 1/12 * fval1
# fval2 <- 1/15 * fval2
# fval3 <- 1/12 * fval3
# fval4 <- 1/10 * fval4
# fval5 <- 1/10 * fval5
# fval6 <- 1/12 * fval6
# fval7 <- 1/8 * fval7
# fval8 <- 1/8 * fval8
# fval9 <- 1/10 * fval9


# basis function check
par(mfrow=c(2,2))
ylim <- quantile( c(fval1, fval2, fval3), c(0,1) )

plot(te, fval1, type='l', ylim=ylim, col=2)
lines(te, fval2, col=3)
lines(te, fval3, col=4)

ylim <- quantile( c(fval4, fval5, fval6), c(0,1) )

plot(te, fval4, type='l', ylim=ylim, col=2)
lines(te, fval5, col=3)
lines(te, fval6, col=4)

ylim <- quantile( c(fval7, fval8, fval9), c(0,1) )

plot(te, fval8, type='l', ylim=ylim, col=2)
lines(te, fval8, col=3)
lines(te, fval9, col=4)



truezs <- c( rep(1,n1), rep(2, n2), rep(3,n3), rep(4,n4), rep(5,n5) )

Ymean1 <- Ymean2 <- Ymean3 <- matrix(NA, length(te), n)
for(i in 1:n){
  
  zi <- truezs[i]
  
  # 1
  # tempu1 <- rnorm(1,    (zi+1),  1)
  # tempu2 <- rnorm(1,    ifelse( zi==2, 1, 0 ),  1)
  # tempu3 <- rnorm(1,    (5-zi),  1)
  # tempu5 <- rnorm(1,    (zi-2),  1)
  # tempu6 <- rnorm(1,    (5-zi),  1)
  # tempu8 <- rnorm(1,    ifelse( zi==3, 1, 0 ),  1)
  # tempu9 <- rnorm(1,    (zi+2),  1)

  # tempu1 <- rnorm(1,    (zi+1),  1.1)
  # tempu2 <- rnorm(1,    ifelse( zi==2, 1, 0 ),  1.1)
  # tempu3 <- rnorm(1,    (5-zi),  1.1)
  # tempu4 <- rnorm(1,    ifelse( zi==3, 1, 0 ),  1.1) 
  # tempu5 <- rnorm(1,    (zi+2),  1.1) 
  # tempu6 <- rnorm(1,    (5-zi),  1.1)
  
  tempu1 <- rnorm(1,    (zi+1),  1.2)
  tempu2 <- rnorm(1,    ifelse( (zi==2||zi==3), -zi+5+1, 0 ) + ifelse( (zi!=2&&zi!=3), zi+1, 0 ),  1.2)
  tempu3 <- rnorm(1,    (6-zi),  1.5)
  tempu4 <- rnorm(1,    (6-zi),  1.2)
  tempu5 <- rnorm(1,    (zi-2),  1.2)
  tempu6 <- rnorm(1,    ifelse( (zi==2||zi==3), -zi+5+1, 0 ) + ifelse( (zi!=2&&zi!=3), zi+1, 0 ),  1.5)
  tempu8 <- rnorm(1,    (zi+2),  1.2)
  tempu9 <- rnorm(1,    (zi-2),  1.2)
  
  tempg11 <- rnorm(1,zi * (-1)^(zi+1),1)
  
  tempg21 <- rnorm(1,zi * (-1)^(zi),1)
  tempg22 <- rnorm(1,(-1)^(zi+1),1)
  
  tempg31 <- rnorm(1,zi * (-1)^(zi),1)
  tempg32 <- rnorm(1,(-1)^(zi),1)
  
  Xdat[,i] <- c(tempu1, tempu2, tempu3, tempu4, tempu5, tempu6, tempu8, tempu9)
  Zdat[,i] <- c( which.max(c(tempg11, 0)), 
                 which.max(c(tempg21, tempg22, 0)),
                 which.max(c(tempg31, tempg32, 0)) )
  
  # # 80
  # Ymean1[,i] <- (zi-1) * (1*te+9) + (tempg11/10)*te + (tempu1+tempg22/10) * fval1 + (tempu2+tempg21/10) * fval2 + tempu3 * fval3
  # Ymean2[,i] <- (zi-2) * (1*te+6) + tempg11/10 + tempu4 * fval4 + (tempu5+tempg31/10) * fval5 + (tempu6+tempg32/10) * fval6 
  # Ymean3[,i] <- (ifelse(zi==3, zi, -zi+3)-2) * (1*te+8) + tempu8 * fval8 + tempu9 * fval9
  
  # # 83
  # Ymean1[,i] <- (zi-1) * (1*te+10) + (tempg11/10)*te + (tempu1+tempg22/10) * fval1 + (tempu2+tempg21/10) * fval2 + tempu3 * fval3
  # Ymean2[,i] <- (zi-2) * (1*te+6) + tempg11/10 + tempu4 * fval4 + (tempu5+tempg31/10) * fval5 + (tempu6+tempg32/10) * fval6 
  # Ymean3[,i] <- (ifelse(zi==3, zi, -zi+3)-2) * (1*te+9) + tempu8 * fval8 + tempu9 * fval9
  
  # 83
  Ymean1[,i] <- (zi-3) * (10*te+8) + (tempg11/10)*te + (tempu1+tempg22/10) * fval1 + (tempu2+tempg21/10) * fval2 + tempu3 * fval3
  Ymean2[,i] <- (ifelse(zi==1, 3, zi-1)-2) * (3*te+9) + tempg11/10 + tempu4 * fval4 + (tempu5+tempg31/10) * fval5 + (tempu6+tempg32/10) * fval6
  Ymean3[,i] <- (ifelse(zi==5, 1, zi+2)-2) * (10*te+5) + tempu8 * fval8 + tempu9 * fval9
}



# for(i in 1:n){ 
#   
#   if(i <= n1){
#     
#     zi <- 1
#     tempu1 <- rnorm(1,    (zi+1),  1.3)
#     tempu2 <- rnorm(1,    (zi+1),  1.3)
#     tempu3 <- rnorm(1,    (5-zi),  1.3)
#     tempu4 <- rnorm(1,    (6-zi),  1.3)
#     tempu5 <- rnorm(1,    (5-zi),  1.3)
#     tempu6 <- rnorm(1,    ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ),  1.3)
#     tempu7 <- rnorm(1,    (zi+2),  1.3)
#     tempu8 <- rnorm(1,    ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ),  1.3)
#     tempu9 <- rnorm(1,    (4-zi),  1.3)
#     
#     
#     
#     # tempu1 <- rnorm(1,    1,  1)
#     # tempu2 <- rnorm(1,    1,  1)
#     # tempu3 <- rnorm(1,    1,  1)
#     # tempu4 <- rnorm(1,    1,  0.5)
#     # tempu5 <- rnorm(1,    0,  1)
#     # tempu6 <- rnorm(1,    0,  1)
#     # 
#     # tempu7 <- rnorm(1,    1,  1)
#     # tempu8 <- rnorm(1,    0,  1)
#     # tempu9 <- rnorm(1,    1,  1)
#     
#     tempg11 <- rnorm(1,1,1)
#     
#     tempg21 <- rnorm(1,-1,1)
#     tempg22 <- rnorm(1,1,1)
#     
#     tempg31 <- rnorm(1,-1,1)
#     tempg32 <- rnorm(1,-1,1)
#     
#     Xdat[,i] <- c(tempu1, tempu2, tempu3, tempu4, tempu5, tempu6, tempu7, tempu8, tempu9)
#     Zdat[,i] <- c( which.max(c(tempg11, 0)), 
#                    which.max(c(tempg21, tempg22, 0)),
#                    which.max(c(tempg31, tempg32, 0)) )
#     
#     Ymean1[,i] <- ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ) + (zi + tempg11/4)*te + (tempu1+tempg22/5) * fval1 + (tempu2+tempg21/5) * fval2 + tempu3 * fval3
#     Ymean2[,i] <- (1.5*zi + tempg11/3) + (tempu4+tempg32/5) * fval4 + (tempu5+tempg31/5) * fval5 + tempu6 * fval6
#     Ymean3[,i] <- (te+1/2)* zi + tempu7 * fval7 + tempu8 * fval8 + tempu9 * fval9
#     
#   }else if(i <= n1+n2){ 
#     
#     zi <- 2
#     tempu1 <- rnorm(1,    (zi+1),  1.3)
#     tempu2 <- rnorm(1,    (zi+1),  1.3)
#     tempu3 <- rnorm(1,    (5-zi),  1.3)
#     tempu4 <- rnorm(1,    (6-zi),  1.3)
#     tempu5 <- rnorm(1,    (5-zi),  1.3)
#     tempu6 <- rnorm(1,    ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ),  1.3)
#     tempu7 <- rnorm(1,    (zi+2),  1.3)
#     tempu8 <- rnorm(1,    ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ),  1.3)
#     tempu9 <- rnorm(1,    (4-zi),  1.3)
#     
#     # tempu1 <- rnorm(1,    2,  1)
#     # tempu2 <- rnorm(1,    2,  1)
#     # tempu3 <- rnorm(1,    2,  1)
#     # tempu4 <- rnorm(1,    0,  0.5)
#     # tempu5 <- rnorm(1,    2,  1)
#     # tempu6 <- rnorm(1,    0,  1)
#     # 
#     # tempu7 <- rnorm(1,    2,  1)
#     # tempu8 <- rnorm(1,    2,  1)
#     # tempu9 <- rnorm(1,    0,  1)
#     
#     tempg11 <- rnorm(1,-2,1)
#     
#     tempg21 <- rnorm(1,2,1)
#     tempg22 <- rnorm(1,-1,1)
#     
#     tempg31 <- rnorm(1,2,1)
#     tempg32 <- rnorm(1,1,1)
#     
#     Xdat[,i] <- c(tempu1, tempu2, tempu3, tempu4, tempu5, tempu6, tempu7, tempu8, tempu9)
#     Zdat[,i] <- c( which.max(c(tempg11, 0)), 
#                    which.max(c(tempg21, tempg22, 0)),
#                    which.max(c(tempg31, tempg32, 0)) )
#     
#     Ymean1[,i] <- ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ) + (zi + tempg11/4)*te + (tempu1+tempg22/5) * fval1 + (tempu2+tempg21/5) * fval2 + tempu3 * fval3
#     Ymean2[,i] <- (1.5*zi + tempg11/3) + (tempu4+tempg32/5) * fval4 + (tempu5+tempg31/5) * fval5 + tempu6 * fval6
#     Ymean3[,i] <- (te+1/2)* zi + tempu7 * fval7 + tempu8 * fval8 + tempu9 * fval9
#     
#   }else{
#     
#     zi <- 3
#     tempu1 <- rnorm(1,    (zi+1),  1.3)
#     tempu2 <- rnorm(1,    (zi+1),  1.3)
#     tempu3 <- rnorm(1,    (5-zi),  1.3)
#     tempu4 <- rnorm(1,    (6-zi),  1.3)
#     tempu5 <- rnorm(1,    (5-zi),  1.3)
#     tempu6 <- rnorm(1,    ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ),  1.3)
#     tempu7 <- rnorm(1,    (zi+2),  1.3)
#     tempu8 <- rnorm(1,    ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ),  1.3)
#     tempu9 <- rnorm(1,    (4-zi),  1.3)
#     
#     # tempu1 <- rnorm(1,    3,  1)
#     # tempu2 <- rnorm(1,    3,  1)
#     # tempu3 <- rnorm(1,    3,  1)
#     # tempu4 <- rnorm(1,    0,  0.5)
#     # tempu5 <- rnorm(1,    0,  1)
#     # tempu6 <- rnorm(1,    3,  1)
#     # 
#     # tempu7 <- rnorm(1,    3,  1)
#     # tempu8 <- rnorm(1,    3,  1)
#     # tempu9 <- rnorm(1,    3,  1)
#     
#     tempg11 <- rnorm(1,3,1)
#     
#     tempg21 <- rnorm(1,-3,1)
#     tempg22 <- rnorm(1,1,1)
#     
#     tempg31 <- rnorm(1,-3,1)
#     tempg32 <- rnorm(1,-1,1)
#     
#     Xdat[,i] <- c(tempu1, tempu2, tempu3, tempu4, tempu5, tempu6, tempu7, tempu8, tempu9)
#     Zdat[,i] <- c( which.max(c(tempg11, 0)), 
#                    which.max(c(tempg21, tempg22, 0)),
#                    which.max(c(tempg31, tempg32, 0)) )
#     
#     Ymean1[,i] <- ifelse( (zi==1), zi+1, 0 ) + ifelse( (zi!=1), 5-zi+1, 0 ) + (zi + tempg11/4)*te + (tempu1+tempg22/5) * fval1 + (tempu2+tempg21/5) * fval2 + tempu3 * fval3
#     Ymean2[,i] <- (1.5*zi + tempg11/3) + (tempu4+tempg32/5) * fval4 + (tempu5+tempg31/5) * fval5 + tempu6 * fval6
#     Ymean3[,i] <- (te+1/2)* zi + tempu7 * fval7 + tempu8 * fval8 + tempu9 * fval9
#   }
# }


Ydat1 <- Ydat11 <- Ymean1 + matrix( rnorm(length(te)*n, 0, 3), length(te), n )
Ydat2 <- Ydat22 <- Ymean2 + matrix( rnorm(length(te)*n, 0, 2), length(te), n )
Ydat3 <- Ydat33 <- Ymean3 + matrix( rnorm(length(te)*n, 0, 4), length(te), n )

Ydat1 <- c(Ydat1)
Ydat2 <- c(Ydat2)
Ydat3 <- c(Ydat3)
Ydat <- list(Ydat1, Ydat2, Ydat3)

tobs1 <- NULL; for(i in 1:n){ tobs1 <- c( tobs1, te ) }
tobs2 <- NULL; for(i in 1:n){ tobs2 <- c( tobs2, te ) }
tobs3 <- NULL; for(i in 1:n){ tobs3 <- c( tobs3, te ) }
tobs <- list(tobs1, tobs2, tobs3)
cc <- 3




######################
# radial basis 
######################
radialBasis <- function( t, Center, nu, p ){ rbasis <- rep(1,p); rbasis[2:p] <- exp( -nu * ( t-Center[-p] )^2); return(rbasis); }

# K2 <- 65
# fb <- create.fourier.basis(c(0,1), K2)
# B2 <- eval.basis(tobs[[2]], fb)

MakeBasis <- function(basis, TOBS, YDAT){
  
  if(basis == 1){
    
    K <- 65
    
    fb <- create.fourier.basis(c(0,1), K)
    B <- eval.basis(TOBS, fb)
    
  }else if(basis == 2){
    
    K <- 40
    
    # nu
    Nu <- 4
    center <- (0:(K-1))/(K-1);
    
    # B
    B <- matrix(0, length(TOBS), K)
    for(i in 1:length(TOBS)){ B[i,] <- radialBasis(TOBS[i], center, Nu, K) }
    
  }else if(basis == 3){
    
    K <- 20
    
    # B
    tempB <- getSplineInfo(TOBS, K)
    K <- ncol(tempB)
    
    Utobs <- sort(unique(TOBS))
    B <- matrix(0, length(TOBS), K)
    for(i in 1:length(TOBS)){ B[i,] <- tempB[which( Utobs %in% TOBS[i] ),] }
    
  }else{
    
    K <- 20
    
    Bbasis <- bs(YDAT, df = K, degree = 3)
    B <- Bbasis
  }
  
  kk <- K
  bb <- B
  
  return( list(kk, bb) )
}

################################################
BasisResult1 <- MakeBasis(3, tobs[[1]], Ydat[[1]]) # 1 2 3 4 
BasisResult2 <- MakeBasis(3, tobs[[2]], Ydat[[2]]) # 1 2 3 4
BasisResult3 <- MakeBasis(3, tobs[[3]], Ydat[[3]]) # 1 2 3 4
################################################

K1 <- BasisResult1[[1]]
B1 <- BasisResult1[[2]]

K2 <- BasisResult2[[1]]
B2 <- BasisResult2[[2]]

K3 <- BasisResult3[[1]]
B3 <- BasisResult3[[2]]

Bcs <- list(B1, B2, B3)
Kcs <- c(K1,K2, K3)
#Be <- getSplineInfo(te, K-4)

# q <- 3
# 
# Zdat <- matrix(NA, q, n)
# 
# Zdat[1,] <- c( rep(1,15), rep(2,15), rep(3,10) )
# Zdat[2,] <- c( rep(1,10), rep(2,10), rep(3,10), rep(4,5), rep(5,5) )
# Zdat[3,] <- rep(1:5, 8)



# par(mfrow=c(2,2), mar=c(5.1, 5.1, 4.1, 2.1))
# 
# YLIM <- quantile(Ydat11, c(0,1))
# plot(te, Ydat11[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
#      main = "(a) bivariate functional data - dimension 1", cex.main=2.5, cex.lab=2 )
# for(i in 1:n1){ lines(te, Ydat11[,i], col=2) }
# for(i in (1:n2)+n1 ){ lines(te, Ydat11[,i],col=3) }
# for(i in (1:n3)+n1+n2 ){ lines(te, Ydat11[,i], col=4) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3"), col=c(2,3,4), lwd=c(2,2,2), cex=2)
# 
# YLIM <- quantile(Ydat22, c(0,1))
# plot(te, Ydat22[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
#      main = "(b) bivariate functional data - dimension 2", cex.main=2.5, cex.lab=2 )
# for(i in 1:n1){ lines(te, Ydat22[,i],col=2) }
# for(i in (1:n2)+n1 ){ lines(te, Ydat22[,i],col=3) }
# for(i in (1:n3)+n1+n2 ){ lines(te, Ydat22[,i], col=4) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3"), col=c(2,3,4), lwd=c(2,2,2), cex=2)
# 
# YLIM <- quantile(Ydat33, c(0,1))
# plot(te, Ydat33[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
#      main = "(b) bivariate functional data - dimension 2", cex.main=2.5, cex.lab=2 )
# for(i in 1:n1){ lines(te, Ydat33[,i],col=2) }
# for(i in (1:n2)+n1 ){ lines(te, Ydat33[,i],col=3) }
# for(i in (1:n3)+n1+n2 ){ lines(te, Ydat33[,i], col=4) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3"), col=c(2,3,4), lwd=c(2,2,2), cex=2)
# 
# 
# 
# 
# figurenum <- c("(a", "(b", "(c", "(d", "(e", "(f", "(g", "(h", "(i", "(j", "(k", "(l", "(m")
# 
# par(mfrow=c(3,3))
# for(i in 1:9){
#   
#   DIM <- i
#   YLIM <- quantile(Xdat[DIM,], c(0,1))
#   
#   plot(1, Xdat[DIM,1], ylim=YLIM, xlim=c(3/4,13/4), xlab='', ylab=paste("u",i,sep=''),
#        main=paste(figurenum[i+4], ") Distribution of ", "U", i, sep=""), xaxt='n', cex.main=1.6, cex.lab=1.3)
#   axis(1, at=c(1,2,3), labels=c("cluster 1", "cluster 2", "cluster 3"))
#   points( 1+rnorm(n1,0,1/15), Xdat[DIM,1:n1] )
#   points( 2+rnorm(n2,0,1/15), Xdat[DIM,1:n2+n1] )
#   points( 3+rnorm(n3,0,1/15), Xdat[DIM,1:n3+n1+n2] )
# }
# 
# for(j in 1:q){
#   
#   x <- table( Zdat[j,], rep(c(1,2,3),c(n1,n2,n3)) )
#   
#   colnames(x) <- c("cluster 1", "cluster 2", "cluster 3")
#   barplot( x, beside=TRUE, ylab='Frequency', col=terrain.colors(length(unique(Zdat[j,]))), legend = rownames(x),
#            args.legend = list(bty = "n", x = "top", ncol = 3), ylim=c(0,24), 
#            main = paste(figurenum[10+j], ") Frequencies of ", "z", j, sep=""), cex.main=1.6, cex.lab=1.3)
# }
# 
# dev.off()


par(mfrow=c(2,2), mar=c(5.1, 5.1, 4.1, 2.1))

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(a) bivariate functional data - dimension 1", cex.main=2.5, cex.lab=2 )
for(i in 1:n1){ lines(te, Ydat11[,i], col=2) }
for(i in (1:n2)+n1 ){ lines(te, Ydat11[,i],col=3) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ydat11[,i], col=4) }
for(i in (1:n4)+n1+n2+n3 ){ lines(te, Ydat11[,i], col=5) }
for(i in (1:n5)+n1+n2+n3+n4 ){ lines(te, Ydat11[,i], col=6) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3"), col=c(2,3,4), lwd=c(2,2,2), cex=2)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(b) bivariate functional data - dimension 2", cex.main=2.5, cex.lab=2 )
for(i in 1:n1){ lines(te, Ydat22[,i],col=2) }
for(i in (1:n2)+n1 ){ lines(te, Ydat22[,i],col=3) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ydat22[,i], col=4) }
for(i in (1:n4)+n1+n2+n3 ){ lines(te, Ydat22[,i], col=5) }
for(i in (1:n5)+n1+n2+n3+n4 ){ lines(te, Ydat22[,i], col=6) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3"), col=c(2,3,4), lwd=c(2,2,2), cex=2)

YLIM <- quantile(Ydat33, c(0,1))
plot(te, Ydat33[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(b) bivariate functional data - dimension 2", cex.main=2.5, cex.lab=2 )
for(i in 1:n1){ lines(te, Ydat33[,i],col=2) }
for(i in (1:n2)+n1 ){ lines(te, Ydat33[,i],col=3) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ydat33[,i], col=4) }
for(i in (1:n4)+n1+n2+n3 ){ lines(te, Ydat33[,i], col=5) }
for(i in (1:n5)+n1+n2+n3+n4 ){ lines(te, Ydat33[,i], col=6) }

