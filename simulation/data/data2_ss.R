library(mvtnorm)
library(fda); library(KFAS); library(MCMCpack); 
#library(stochvol);
library(truncnorm) ; 
library(dlm)
library(MixMatrix)
library(CVTuningCov)

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
nts <- rep(nn, TT)
te <- seq(0,1,length.out = ss) #seq(1,21,0.05)

# observation numbers
nobs1 <- nobs2 <- rep(length(te),sum(nts)) #sample(40:60, n, replace=TRUE)
nobs <- list(nobs1, nobs2)





rho <- 0.8
VV <- AR1(TT,rho=rho)

trueprobs <- rmatrixnorm(n = 1, mean=matrix(3,5,TT), U=0.5*diag(5), V=VV)
trueprobs <- sapply(1:TT, function(x){ rdirichlet(1, trueprobs[,x]) })
trueprobs

truezs <- NULL
for(tt in 1:TT){ truezs <- c(truezs, apply(rmultinom(nn,1,trueprobs[,tt]), 2, function(x){ which( x==1 ) })) }


p <- 6
Xdat <- matrix(NA, p, sum(nts))


# fval1 <- 4*log( te + 1/10 ) + 5
# fval2 <- 10*cos(5*te)+2
# fval3 <- exp(2*te)
# fval4 <- tan(te+1/4)
# fval5 <- 10*sin( 3*te + 1/2 ) - 5
# fval6 <- 5*(te-1)^2 - 10
# 
# fval1 <- 1/5 * fval1
# fval2 <- 1/5 * fval2
# fval3 <- 1/5 * fval3
# fval4 <- 1/5 * fval4
# fval5 <- 1/5 * fval5
# fval6 <- 1/5 * fval6


fval1 <- 2*log( te + 1/10 )+2
fval2 <- 5/2*cos(5*te)+2 - 3/2
fval3 <- exp(te)
fval4 <- tan(te+1/3)
fval5 <- 5*sin( 3*te + 1/2 ) - 2
fval6 <- te^3 - 3*te^2


# basis function check
par(mfrow=c(2,2))
ylim <- quantile( c(fval1, fval2, fval3), c(0,1) )

plot(te, fval1, type='l', ylim=ylim, col=1)
lines(te, fval2, col=2)
lines(te, fval3, col=3)

ylim <- quantile( c(fval4, fval5, fval6), c(0,1) )

plot(te, fval4, type='l', ylim=ylim, col=1)
lines(te, fval5, col=2)
lines(te, fval6, col=3)




Ymean1 <- Ymean2 <- matrix(NA, length(te), sum(nts))
for(i in 1:sum(nts)){ 
  
  zi <- truezs[i]
  
  tempu1 <- rnorm(1,    (zi+1),  0.7)
  tempu2 <- rnorm(1,    ifelse( (zi==2&&zi==3), -zi+5+1, 0 ) + ifelse( (zi!=2&&zi!=3), zi+1, 0 ),  1.2)
  tempu3 <- rnorm(1,    (6-zi),  0.7)
  tempu4 <- rnorm(1,    (6-zi),  0.7)
  tempu5 <- rnorm(1,    (zi-2),  0.7)
  tempu6 <- rnorm(1,    ifelse( (zi==2&&zi==3), -zi+5+1, 0 ) + ifelse( (zi!=2&&zi!=3), zi+1, 0 ),  1.2)
  
  Xdat[,i] <- c(tempu1, tempu2, tempu3, tempu4, tempu5, tempu6)
  
  Ymean1[,i] <- (3*te+5)*ifelse( (zi!=2&&zi!=4), zi, -zi ) + tempu1 * fval1 + tempu2 * fval2 + tempu3 * fval3
  Ymean2[,i] <- ((10)*te+7)*zi + tempu4 * fval4 + tempu5 * fval5 + tempu6 * fval6
}


Ydat1 <- Ydat11 <- Ymean1 + matrix( rnorm(length(te)*sum(nts), 0, 4), length(te), sum(nts) )
Ydat2 <- Ydat22 <- Ymean2 + matrix( rnorm(length(te)*sum(nts), 0, 4), length(te), sum(nts) )

Ydat1 <- c(Ydat1)
Ydat2 <- c(Ydat2)
Ydat <- list(Ydat1, Ydat2)

tobs1 <- NULL; for(i in 1:sum(nts)){ tobs1 <- c( tobs1, te ) }
tobs2 <- NULL; for(i in 1:sum(nts)){ tobs2 <- c( tobs2, te ) }
tobs <- list(tobs1, tobs2)
cc <- 2







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
################################################

K1 <- BasisResult1[[1]]
B1 <- BasisResult1[[2]]

K2 <- BasisResult2[[1]]
B2 <- BasisResult2[[2]]

Bcs <- list(B1, B2)

#Be <- getSplineInfo(te, K-4)
n <- sum(nts)

Kcs <- c(K1,K2)
















par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))

tt <- 1
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(a) bivariate functional data - dimension 1 at t=1", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(b) bivariate functional data - dimension 2 at t=1", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)



tt <- 2
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(c) bivariate functional data - dimension 1 at t=2", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(d) bivariate functional data - dimension 2 at t=2", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)



tt <- 3
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(e) bivariate functional data - dimension 1 at t=3", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(f) bivariate functional data - dimension 2 at t=3", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)



tt <- 4
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(g) bivariate functional data - dimension 1 at t=4", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(h) bivariate functional data - dimension 2 at t=4", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)








par(mfrow=c(2,3))

i <- 1
DIM <- i
temp <- Xdat[DIM,]
temp2 <- truezs
temp2[truezs == 1] <- "Cluster 1"
temp2[truezs == 2] <- "Cluster 2"
temp2[truezs == 3] <- "Cluster 3"
temp2[truezs == 4] <- "Cluster 4"
temp2[truezs == 5] <- "Cluster 5"
temptemp <- data.frame(x=temp, y=as.factor(temp2))
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[1]), xlab="", main=expression("(i) Distribution of U"[1]), cex.main=1.6, cex.lab=1.3)

i <- 2
DIM <- i
temp <- Xdat[DIM,]
temp2 <- truezs
temp2[truezs == 1] <- "Cluster 1"
temp2[truezs == 2] <- "Cluster 2"
temp2[truezs == 3] <- "Cluster 3"
temp2[truezs == 4] <- "Cluster 4"
temp2[truezs == 5] <- "Cluster 5"
temptemp <- data.frame(x=temp, y=as.factor(temp2))
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[2]), xlab="", main=expression("(j) Distribution of U"[2]), cex.main=1.6, cex.lab=1.3)

i <- 3
DIM <- i
temp <- Xdat[DIM,]
temp2 <- truezs
temp2[truezs == 1] <- "Cluster 1"
temp2[truezs == 2] <- "Cluster 2"
temp2[truezs == 3] <- "Cluster 3"
temp2[truezs == 4] <- "Cluster 4"
temp2[truezs == 5] <- "Cluster 5"
temptemp <- data.frame(x=temp, y=as.factor(temp2))
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[3]), xlab="", main=expression("(k) Distribution of U"[3]), cex.main=1.6, cex.lab=1.3)

i <- 4
DIM <- i
temp <- Xdat[DIM,]
temp2 <- truezs
temp2[truezs == 1] <- "Cluster 1"
temp2[truezs == 2] <- "Cluster 2"
temp2[truezs == 3] <- "Cluster 3"
temp2[truezs == 4] <- "Cluster 4"
temp2[truezs == 5] <- "Cluster 5"
temptemp <- data.frame(x=temp, y=as.factor(temp2))
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[4]), xlab="", main=expression("(i) Distribution of U"[4]), cex.main=1.6, cex.lab=1.3)

i <- 5
DIM <- i
temp <- Xdat[DIM,]
temp2 <- truezs
temp2[truezs == 1] <- "Cluster 1"
temp2[truezs == 2] <- "Cluster 2"
temp2[truezs == 3] <- "Cluster 3"
temp2[truezs == 4] <- "Cluster 4"
temp2[truezs == 5] <- "Cluster 5"
temptemp <- data.frame(x=temp, y=as.factor(temp2))
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[5]), xlab="", main=expression("(j) Distribution of U"[5]), cex.main=1.6, cex.lab=1.3)

i <- 6
DIM <- i
temp <- Xdat[DIM,]
temp2 <- truezs
temp2[truezs == 1] <- "Cluster 1"
temp2[truezs == 2] <- "Cluster 2"
temp2[truezs == 3] <- "Cluster 3"
temp2[truezs == 4] <- "Cluster 4"
temp2[truezs == 5] <- "Cluster 5"
temptemp <- data.frame(x=temp, y=as.factor(temp2))
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[6]), xlab="", main=expression("(k) Distribution of U"[6]), cex.main=1.6, cex.lab=1.3)

dev.off()










par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))

tt <- 1
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(a) bivariate functional data - dimension 1 at t=1", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(b) bivariate functional data - dimension 2 at t=1", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
# legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)



tt <- 2
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(c) bivariate functional data - dimension 1 at t=2", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(d) bivariate functional data - dimension 2 at t=2", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
# legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)



tt <- 3
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(e) bivariate functional data - dimension 1 at t=3", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(f) bivariate functional data - dimension 2 at t=3", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
# legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)



tt <- 4
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(g) bivariate functional data - dimension 1 at t=4", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
# legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(h) bivariate functional data - dimension 2 at t=4", cex.main=3, cex.lab=2)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
# legend("bottomleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.8)



