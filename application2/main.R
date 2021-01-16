####################################################################################
# Nonparametric Bayesian latent factor model for multivariate functional data
####################################################################################
# by Daewon Yang
####################################################################################
# 2020.09.13
####################################################################################
# Japanese suicide seasonality
####################################################################################
rm(list=ls())


##########################
# load library
##########################
library(mvtnorm)
library(fda)
library(MCMCpack)
library(Rcpp)
library(RcppArmadillo)
library(maps)
library(mondate)
library(mgcv)
library(lubridate)
library(mixmeta)


##########################
# data
##########################
load("Japan8612.RData")

citynames <- city$cityname
n <- length(citynames)

dt <- data[[1]]
days <- format(dt$date[year(dt$date) == year(dt$date)[1]], format="%m-%d")
years <- sort(unique(year(dt$date)))
oneyear_tobs <- (0:364)/364

tobs1 <- tobs2 <- nobs1 <- nobs2 <- Ydat1 <- Ydat2 <- NULL

for(i in 1:n){
  
  cat(i,"\n")
  
  dt <- data[[i]]
  
  onecity_tmp_mat <- matrix(NA, length(years), 365)
  onecity_sun_mat <- matrix(NA, length(years), 365)
  
  for(y in years){
    
    tempdate <- dt$date[year(dt$date) == y]
    tempdate <- format(tempdate, format="%m-%d")
    
    # delete 02-29
    temptmean <- dt$tmean[year(dt$date) == y]
    temptmean <- temptmean[tempdate %in% days]
    
    tempsun <- dt$sun[year(dt$date) == y]
    tempsun <- tempsun[tempdate %in% days]
    
    if( sum(tempdate %in% "02-29") ){ 
      
      tempdate <- tempdate[!(tempdate %in% "02-29")] 
      temptmean <- temptmean[!(tempdate %in% "02-29")]
      tempsun <- tempsun[!(tempdate %in% "02-29")]
    }
    
    # delete NA
    naind <- which(is.na(temptmean))
    if(length(naind)>0){
      
      tempdate <- tempdate[!naind]
      temptmean <- temptmean[!naind]
      tempsun <- tempsun[!naind]
    }
    
    onecity_tmp_mat[which(years==y),days %in% tempdate] <- temptmean
    onecity_sun_mat[which(years==y),days %in% tempdate] <- tempsun
  }
  
  onecity_tmp <- apply(onecity_tmp_mat, 2, function(x){ mean(x,na.rm=T) })
  onecity_sun <- apply(onecity_sun_mat, 2, function(x){ mean(x,na.rm=T) })
  dayind <- 1:365
  
  # delete NA
  naind <- which(is.na(onecity_tmp) | is.na(onecity_sun))
  if(length(naind)>0){
    
    onecity_tmp <- onecity_tmp[!naind]
    onecity_sun <- onecity_sun[!naind]
    dayind <- dayind[!naind]
  }
  
  tobs1 <- c(tobs1, oneyear_tobs[dayind])
  nobs1 <- c(nobs1, length(onecity_tmp))
  Ydat1 <- c(Ydat1, onecity_tmp)
  
  tobs2 <- c(tobs2, oneyear_tobs[dayind])
  nobs2 <- c(nobs2, length(onecity_sun))
  Ydat2 <- c(Ydat2, onecity_sun)
}

Ydat <- cbind(Ydat1, Ydat2)
tobs <- cbind(tobs1, tobs2)
nobs <- cbind(nobs1, nobs2)

# load("Japan7215.RData")


city_longi <- city$long.x
city_lati <- city$lat.x
spartial_order <- c(1, 2, 3, 4, 5, 35, 7, 8, 36, 10, 11, 12, 13, 46, 15, 31, 17, 40, 19, 
                    20, 21, 22, 23, 24, 25, 26, 27, 28, 42, 30, 16, 32, 33, 34, 6, 9, 37, 
                    38, 14, 18, 41, 29, 39, 44, 45, 43, 47)

city_longi <- city_longi[spartial_order]
city_lati <- city_lati[spartial_order]



##########################
# exploratory analysis
##########################
monthnum <- c(16, 47, 75, 106, 136, 167, 197,  228, 259, 289, 320, 350)
monthnames <- c("Jan.", "Feb.", "Mar.", "Apr.", "May", "June", "July", "Aug.", "Sep.", "Oct.", "Nov.", "Dec.")


# temperature
tempind1 <- cumsum(c(0,nobs1))

plot( 1:365, Ydat1[(tempind1[1]+1):tempind1[2]], type='l', ylim=range(Ydat1), main="Temperature in Japan", xlab='Month', ylab='Temperature', xaxt='n' )
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]] ) }


# sunshine
tempind2 <- cumsum(c(0,nobs2))

plot( 1:365, Ydat2[(tempind2[1]+1):tempind2[2]], type='l', ylim=range(Ydat2), main="Sunshine in Japan", xlab='Month', ylab='Sunshine', xaxt='n' )
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:365, Ydat2[(tempind2[i]+1):tempind2[i+1]] ) }


### plot 1
library(NipponMap)
aaa <- rep(0,n)
tempcityind <- c(1,21,47)
tempcitynames <- citynames[tempcityind]
aaa[tempcityind] <- c(2,3,4)
a2 <- JapanPrefMap(aaa, border = 1, axes = TRUE)

text( a2[tempcityind[1],1]+3, a2[tempcityind[1],2]-1, tempcitynames[1], cex=1, col=2 )
segments( a2[tempcityind[2],1], a2[tempcityind[2],2], a2[tempcityind[2],1]+0.5, a2[tempcityind[2],2]-2, col=3  )
text( a2[tempcityind[2],1]+0.6, a2[tempcityind[2],2]-2.3, tempcitynames[2], cex=1, col=3 )
text( 135, 41.3, tempcitynames[3], cex=1, col=4 )
title(main="(a)", cex.main=1.5)

# plot 2
plot( 1:365, Ydat1[(tempind1[1]+1):tempind1[2]], type='l', ylim=range(Ydat1), main="(b) Temperature in Japan", xlab='Month', ylab='Temperature (¡ÆC)', xaxt='n', cex.main=1.5, cex.lab=1.3 )
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)

i <- tempcityind[1]
lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]], col=2, lwd=1.5 )
i <- tempcityind[2]
lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]], col=3, lwd=1.5 )
i <- tempcityind[3]
lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]], col=4, lwd=1.5 )
legend("topleft", tempcitynames, lwd=c(2,2,2), col=c(2,3,4), cex=1.4)

# plot 3
plot( 1:365, Ydat2[(tempind2[1]+1):tempind2[2]], type='l', ylim=c(1.5,9.5), main="(c) Sunshine in Japan", xlab='Month', ylab='Sunshine (hour)', xaxt='n', cex.main=1.5, cex.lab=1.3 )
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)

i <- tempcityind[1]
lines( 1:365, Ydat2[(tempind1[i]+1):tempind1[i+1]], col=2, lwd=1.5 )
i <- tempcityind[2]
lines( 1:365, Ydat2[(tempind1[i]+1):tempind1[i+1]], col=3, lwd=1.5 )
i <- tempcityind[3]
lines( 1:365, Ydat2[(tempind1[i]+1):tempind1[i+1]], col=4, lwd=1.5 )
legend("topright", tempcitynames, lwd=c(2,2,2), col=c(2,3,4), cex=1.4)


##########################
# Two-stage modeling
##########################
source("cyclic.R")

# degree of freedom 
dfval <- 4 
knots7 <- seq(1+52/8, 53-52/8, length=7)
knots4 <- knots7[c(1,3,5,7)]


###### first stage #######
NumZero<-matrix(NA,n,1)

data_week <- list()
VCOV <- list()
COEF <- NULL

for(i in 1:n){
  
  dt <- data[[i]]
  
  ## Numbering week 1. 01-01 = start of week-of-year
  dt$year <- as.factor(format(dt$date, format="%Y"))
  dt$week <- as.factor(week(dt$date))
  dt$yrwk <- as.factor(dt$year:dt$week)
  
  # weekly data
  dt_week<-data.frame(year=as.factor(rep(unique(dt$year),each=53)),
                      week=unique(week(dt$date)),
                      stot=tapply(dt$stot,dt$yrwk,sum))
  
  # correction
  tempsuic2 <- as.numeric(tapply(dt$date, as.factor(as.factor(dt$yr):as.factor(week(dt$date))), length))
  wk53ind <- which(tempsuic2 != 7)
  dt_week$stot[wk53ind] <- sapply( wk53ind, function(i){ floor(sum(dt_week$stot[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]],na.rm=T) * (7/sum(tempsuic2[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]]))) } )
  dt_week$week1 <- as.numeric(dt_week$week == 1)
  
  # number of zero stratum
  StotStra<-tapply(dt_week$stot,dt_week$year,sum)
  drop.zero<-names(which(StotStra==0))
  
  city_ind <- i
  NumZero[city_ind]<-length(which(dt_week$year %in% drop.zero))
  
  # remove stratum having zero suicide count
  if(NumZero[city_ind]>0){   dt_week<-dt_week[-which(dt_week$year %in% drop.zero),]}
  
  # Weekly data set without zero stratum
  data_week[[city_ind]]<-dt_week
  
  ### GLM
  # seasonality basis (cyclic spline)
  temp_spline.season <- cyclic(dt_week$week,df=4, knots=knots4)
  
  model <- glm(stot ~ temp_spline.season + year, data=dt_week, family=quasipoisson)
  # COEF[[city_ind]] <- model$coefficients[2:(dfval+1)]
  COEF <- rbind(COEF, model$coefficients[2:(dfval+1)])
  VCOV[[city_ind]] <- vcov(model)[2:(dfval+1),2:(dfval+1)]
}

###### second stage #######
metaRESULT <- mixmeta(COEF, S=VCOV)
BLUPresult <- blup(metaRESULT)


basisval <- cyclic(1:53, df=4, knots=knots4)
logRRval <- basisval %*% t(BLUPresult)
logRRval <- logRRval - t(matrix( logRRval[1,], n, 53 ))
RRval <- exp(logRRval)

plot(1:53, RRval[,1], type='l', ylim=range(RRval), xlab="Month", ylab="RR", xaxt='n')
axis(1, c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5, 
          3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), monthnames, labels = F)
text( c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5, 
        3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:53, RRval[,i] ) }





##########################
# main - clustering
##########################
sourceCpp("fastNEW2.cpp")

#### basis function
radialBasis <- function( t, Center, nu, p ){ rbasis <- rep(1,p); rbasis[2:p] <- exp( -nu * ( t-Center[-p] )^2); return(rbasis); }

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


MakeBasis <- function(K, basis, TOBS, YDAT){
  
  if(basis == 1){
    
    # K <- 65
    
    fb <- create.fourier.basis(c(0,1), K)
    B <- eval.basis(TOBS, fb)
    
  }else if(basis == 2){
    
    # K <- 40
    
    # nu
    Nu <- 4
    center <- (0:(K-1))/(K-1);
    
    # B
    B <- matrix(0, length(TOBS), K)
    for(i in 1:length(TOBS)){ B[i,] <- radialBasis(TOBS[i], center, Nu, K) }
    
  }else if(basis == 3){
    
    # K <- 20
    
    # B
    tempB <- getSplineInfo(TOBS, K)
    K <- ncol(tempB)
    
    Utobs <- sort(unique(TOBS))
    B <- matrix(0, length(TOBS), K)
    for(i in 1:length(TOBS)){ B[i,] <- tempB[which( Utobs %in% TOBS[i] ),] }
    
  }else{
    
    # K <- 20
    
    Bbasis <- bs(YDAT, df = K, degree = 3)
    B <- Bbasis
  }
  
  kk <- K
  bb <- B
  
  return( list(kk, bb) )
}

niter <- 15000; 
burn <- 10000
thin <- 1
burnthin <- seq(burn+1, niter, thin)
# num <- 20


basislist <- c(1,2,3,4)
Kvals <- c(65,40,20,20)

b <- 3
basis <- basislist[b]
K <- Kvals[b]

# basis
Kcs <- NULL
Bcs <- list()

for(i in 1:ncol(Ydat)){
  
  tempBasisResult <- MakeBasis(K, basis, tobs[,i], Ydat[,i])
  
  tempK <- tempBasisResult[[1]]
  tempB <- tempBasisResult[[2]]
  
  Kcs <- c(Kcs, tempK)
  Bcs[[i]] <- tempB
}

# data for suicide seasonality
H <- dfval
hatPSIis <- t(COEF)
Svaris <- NULL
for(i in 1:n){ Svaris <- cbind(Svaris, VCOV[[i]]) }

# without covariate
containCOV <- 1
containCAT <- 0
containMULT <- 1

Xdat <- rbind(city_longi, city_lati)

set.seed(5) #  1 -> 120 5 100 2 4 17 19
             # 0.5 -> 30
source("preprocess2.R")

LL_result <- matrix(0, length(burnthin), n)
PS_Zmu <- list()
PS_LFLM <- list()

for(wert in 1:niter){
  
  cat(wert, "th iteration \n")
  
  ### part 1
  tempzv = getZandV( n, M, alpha, beta, tempKHpq, ppi, llamb, ZETA, a_u, b_u, a_w, b_w, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu )
  ZZ = tempzv[1:tempKHpq,]
  VV = tempzv[1:tempKHpq+tempKHpq,]
  Eta = t( tempzv[1:n+2*tempKHpq,] )
  ums = c( tempzv[2*tempKHpq+n+1,] )
  mm = c( tempzv[2*tempKHpq+n+2,] )
  wms = c( tempzv[2*tempKHpq+n+3,] )
  Zmu = t( tempzv[1:M+2*tempKHpq+n+3,] )
  
  LFLM = ZZ * VV
  r = nrow(Eta)
  
  cat("r = ", r, "\n")
  
  ums = getUMS(r, ZZ, VV, a_u, b_u);
  Eta = getEta(n, r, LFLM, Sigmas, ZETA, LL, Zmu );
  
  alpha = getAlpha(ZZ, containBETA, tempKHpq, beta, a_alpha, b_alpha, Hp);
  beta = getBeta(ZZ, containBETA, tempKHpq, a_beta, b_beta, beta, alpha);
  
  ### part 2
  Zmu = getZmu(LL, Eta, r, wms, mm, M);
  PP = getPIS(M, eps, LL, nu);
  LL = getZis( n, M, r, PP, Eta, Zmu );
  mm = getmm(M, r, wms, Zmu);    
  wms = getwms(M, Zmu, a_w, b_w, mm, r);
  nu = getdelta( M, PP, a_nu, b_nu );
  
  ### part 3
  Sigmas = getSigmas( tempKHpq, a_sigma, b_sigma, Eta, ZZ, VV, ZETA, n );
  
  if( containCAT == 1 ){ gis = getgis( tempKHp, qind1, qind, qdims, Zdatstar, gis, n, LFLM, Eta, Sigmas, qq ); }
  if( cc>1 ){
    
    tempTHETA = getTHETA( cc, kkind, n, nobs, B, Sigmas, KKs, psi, LFLM, Eta, Ydat );
    psi = getpsi( n, a_psi, b_psi, cc, nobs, KKs, tempTHETA, Ydat, B );
    
    tempPSIis <- getPSIis( n, Sigmas, tempK, dfval, LFLM, Eta, hatPSIis, Svaris )
    
    if(containCOV == 1){
      
      ZETA = rbind(tempTHETA, tempPSIis)
      ZETA = rbind(ZETA, Xdat)
      
    }else{
      
      ZETA = rbind(tempTHETA, tempPSIis)
    }
    
  }else{
    
    tempTHETA = getTHETA0( n, nobs, B, Sigmas, tempK, psi0, LFLM, Eta, Ydat );
    psi0 = getpsi0( n, a_psi, b_psi, nobs, tempTHETA, Ydat, B );
    
    tempPSIis <- getPSIis( n, Sigmas, tempK, dfval, LFLM, Eta, hatPSIis, Svaris )
    
    if(containCOV == 1){
      
      ZETA = rbind(tempTHETA, tempPSIis)
      ZETA = rbind(ZETA, Xdat)
      
    }else{
      
      ZETA = rbind(tempTHETA, tempPSIis)
    }
  }
  
  if( containCAT == 1 ){ ZETA = rbind( ZETA, gis ); }
  
  
  if( wert %in% burnthin ){
    
    LL_result[which(burnthin==wert),] = LL
    PS_Zmu[[which(burnthin %in% wert)]] = Zmu
    PS_LFLM[[which(burnthin %in% wert)]] = LFLM
  }
}

clust_result <- Dahlclust(LL_result)
clusters <- clust_result[,1]
tempresult2 <- clusters



library(plot.matrix)
order(clusters)
hatmat <- clust_result[,-1]
plot(hatmat[order(clusters),order(clusters)])


# temperature
tempind1 <- cumsum(c(0,nobs1))

plot( 1:365, Ydat1[(tempind1[1]+1):tempind1[2]], type='l', ylim=range(Ydat1), main="Temperature in Japan", xlab='Month', ylab='Temperature', xaxt='n' )
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]], col=tempresult2[i] ) }

# sunshine
tempind2 <- cumsum(c(0,nobs2))

plot( 1:365, Ydat2[(tempind2[1]+1):tempind2[2]], type='l', ylim=range(Ydat2), main="Sunshine in Japan", xlab='Month', ylab='Sunshine', xaxt='n' )
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:365, Ydat2[(tempind2[i]+1):tempind2[i+1]], col=tempresult2[i] ) }


basisval <- cyclic(1:53, df=4, knots=knots4)
logRRval <- basisval %*% t(BLUPresult)
logRRval <- logRRval - t(matrix( logRRval[1,], n, 53 ))
RRval <- exp(logRRval)

plot(1:53, RRval[,1], type='l', ylim=range(RRval), xlab="Month", ylab="RR", xaxt='n')
axis(1, c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5,
          3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), monthnames, labels = F)
text( c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5,
        3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:53, RRval[,i], col=tempresult2[i] ) }





table(clusters)






#####################################
# plot
#####################################

# cluster specific mean
nind <- c(0,cumsum(nobs))

DAT <- Ydat
DIM <- 1
tempDAT <- DAT[,DIM]
DIM2 <- 2
tempDAT2 <- DAT[,DIM2]
YLIM <- quantile( tempDAT, c(0,1) )
YLIM2 <- quantile( tempDAT2, c(0,1) )

library(abind)
library(label.switching)
tempKval <- max(LL_result)
run <- ecr(zpivot = clusters, z = LL_result, K = tempKval)

THETAmean_star <- array(NA, c(nrow(LL_result), tempK, tempKval) )
THETAmean_star2 <- array(NA, c(nrow(LL_result), dfval, tempKval) )
for(i in 1:nrow(LL_result)){
  
  perm <- run$permutations[i,]
  tempLL <- LL_result[i,]
  tempZmu <- PS_Zmu[[i]]
  tempZmu2 <- tempZmu[,1:tempKval]
  tempLFLM <- PS_LFLM[[i]]
  
  THETAmean_star[i,,] <- (tempLFLM %*% (tempZmu2[,perm]))[1:tempK,]
  THETAmean_star2[i,,] <- (tempLFLM %*% (tempZmu2[,perm]))[1:dfval+tempK,]
}

clusters
THETAmean_star <- apply(THETAmean_star, c(2,3), mean)
THETAmean_star2 <- apply(THETAmean_star2, c(2,3), mean)

ii <- 1
msc1 <- Bcs[[1]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_star[1:Kcs[1],]
msc2 <- Bcs[[2]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_star[1:Kcs[2]+Kcs[1],]

basisval <- cyclic(1:53, df=4, knots=knots4)
msc3 <- basisval %*% THETAmean_star2


clusters2 <- clusters
clusters2[clusters==2] <- 4
clusters2[clusters==4] <- 3
clusters2[clusters==3] <- 2

# clusters2[clusters==2] <- 4
# clusters2[clusters==3] <- 2
# clusters2[clusters==4] <- 1
# clusters2[clusters==3] <- 4

# 1 2 3 5
# 6
# 47


monthnum <- c(16, 47, 75, 106, 136, 167, 197,  228, 259, 289, 320, 350)
monthnum2 <- c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5, 3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 )
monthnames <- c("Jan.", "Feb.", "Mar.", "Apr.", "May", "June", "July", "Aug.", "Sep.", "Oct.", "Nov.", "Dec.")

plot(msc1[,1], type='l', ylim=c(-2,30), main="(a) Temperature", ylab="Temperature (¡ÆC)", xlab='Month', xaxt='n', cex.main=1.5, cex.lab=1.3)
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }


plot(msc2[,1], type='l', ylim=c(0,10), main="(b) Sunshine", ylab="Sunshine", xlab='Month', xaxt='n', cex.main=1.5, cex.lab=1.3)
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }

plot(msc3[,1], type='l', ylim=c(-0.2,0.01), main="(c) logRR", ylab="logRR", xlab='Month', xaxt='n')
axis(1, monthnum2, monthnames, labels = F)
text(monthnum2, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc3[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }

msc33 <- msc3 - t(matrix(msc3[1,], tempKval, 53))
plot(msc33[,1], type='l', ylim=c(-0.01,0.2), main="(c) logRR", ylab="logRR", xlab='Month', xaxt='n')
axis(1, monthnum2, monthnames, labels = F)
text(monthnum2, par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc33[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }

msc33 <- msc3 - t(matrix(msc3[1,], tempKval, 53))
plot(exp(msc33[,1]), type='l', ylim=c(1,1.32), main="(c) RR", ylab="RR", xlab='Month', xaxt='n', cex.main=1.5, cex.lab=1.3)
axis(1, monthnum2, monthnames, labels = F)
text(monthnum2, par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(exp(msc33[,clusters[which(clusters2 == i)[1]]]), col=i, lwd=3) }


# map
library(NipponMap)

par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
JapanPrefMap(clusters2, border = gray(.8), axes = TRUE)
title(main="(c)", cex.main=1.5)
legend(150, 45.56, legend = c("cluster 1", "cluster 2", "cluster 3"), col=c(1,2,4), pch = c(3,3,3), cex=1.5)




par(mfrow=c(2,2))

plot(msc1[,2], type='l', ylim=c(-2,30), main="(a) Temperature", ylab="Temperature (¡ÆC)", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }


plot(msc2[,2], type='l', ylim=c(2,7), main="(b) Sunshine", ylab="Sunshine (hour)", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }

msc33 <- msc3 - t(matrix(msc3[1,], tempKval, 53))
plot(exp(msc33[,2]), type='l', ylim=c(1,1.35), main="(c) RR", ylab="RR", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum2, monthnames, labels = F)
text(monthnum2, par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(exp(msc33[,clusters[which(clusters2 == i)[1]]]), col=i, lwd=3) }

par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
JapanPrefMap(clusters2, border = gray(.8), axes = TRUE)
title(main="(d)", cex.main=1.7)
legend(153.1, 45.56, legend = c("cluster 1", "cluster 2", "cluster 3"), col=c(2,3,4), pch = c(3,3,3), cex=1.6)

dev.off()




par(mfrow=c(3,3))

msc33 <- msc3 - t(matrix(msc3[1,], tempKval, 53))
plot(exp(msc33[,1]), type='l', ylim=c(1,1.32), main="(c) RR", ylab="RR", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum2, monthnames, labels = F)
text(monthnum2, par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(exp(msc33[,clusters[which(clusters2 == i)[1]]]), col=i, lwd=3) }

plot(msc3[,1], type='l', ylim=range(msc3), main="(c) RR", ylab="RR", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum2, monthnames, labels = F)
text(monthnum2, par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc3[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }


basisval <- cyclic(1:53, df=4, knots=knots4)
logRRval <- basisval %*% t(COEF)

plot(1:53, logRRval[,1], type='l', ylim=range(logRRval), xlab="Month", ylab="RR", xaxt='n')
axis(1, c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5,
          3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), monthnames, labels = F)
text( c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5,
        3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:53, logRRval[,i], col=tempresult2[i] ) }

logRRval <- logRRval - t(matrix( logRRval[1,], n, 53 ))
RRval <- exp(logRRval)

plot(1:53, RRval[,1], type='l', ylim=range(RRval), xlab="Month", ylab="RR", xaxt='n')
axis(1, c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5,
          3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), monthnames, labels = F)
text( c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5,
        3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
for(i in 1:n){ lines( 1:53, RRval[,i], col=tempresult2[i] ) }

hatmat <- clust_result[,-1]
plot(hatmat[order(clusters),order(clusters)])

JapanPrefMap(clusters2, border = gray(.8), axes = TRUE)
title(main="(d)", cex.main=1.7)
legend(152.5, 45.56, legend = c("cluster 1", "cluster 2", "cluster 3"), col=c(1,2,4), pch = c(3,3,3), cex=1.7)

plot(msc1[,1], type='l', ylim=c(-2,30), main="(a) Temperature", ylab="Temperature (¡ÆC)", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }


plot(msc2[,1], type='l', ylim=c(0,10), main="(b) Sunshine", ylab="Sunshine", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }

msc33 <- msc3 - t(matrix(msc3[1,], tempKval, 53))
plot(exp(msc33[,1]), type='l', ylim=c(1,1.32), main="(c) RR", ylab="RR", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
axis(1, monthnum2, monthnames, labels = F)
text(monthnum2, par("usr")[3] - 0.01, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(exp(msc33[,clusters[which(clusters2 == i)[1]]]), col=i, lwd=3) }




















# 
# 
# 
# ### plot 1
# library(NipponMap)
# aaa <- rep(0,n)
# tempcityind <- c(1,21,47)
# tempcitynames <- citynames[tempcityind]
# aaa[tempcityind] <- c(2,3,4)
# a2 <- JapanPrefMap(aaa, border = 1, axes = TRUE)
# 
# text( a2[tempcityind[1],1]+3, a2[tempcityind[1],2]-1, tempcitynames[1], cex=1, col=2 )
# segments( a2[tempcityind[2],1], a2[tempcityind[2],2], a2[tempcityind[2],1]+0.5, a2[tempcityind[2],2]-2, col=3  )
# text( a2[tempcityind[2],1]+0.6, a2[tempcityind[2],2]-2.3, tempcitynames[2], cex=1, col=3 )
# text( 135, 41.3, tempcitynames[3], cex=1, col=4 )
# title(main="(a)", cex.main=1.5)
# 
# # plot 2
# plot( 1:365, Ydat1[(tempind1[1]+1):tempind1[2]], type='l', ylim=range(Ydat1), main="(c) Temperature in Japan", xlab='Month', ylab='Temperature (¡ÆC)', xaxt='n', cex.main=1.5, cex.lab=1.3 )
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
# 
# i <- tempcityind[1]
# lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]], col=2, lwd=1.5 )
# i <- tempcityind[2]
# lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]], col=3, lwd=1.5 )
# i <- tempcityind[3]
# lines( 1:365, Ydat1[(tempind1[i]+1):tempind1[i+1]], col=4, lwd=1.5 )
# legend("topleft", tempcitynames, lwd=c(2,2,2), col=c(2,3,4), cex=1.4)
# 
# # plot 3
# plot( 1:365, Ydat2[(tempind2[1]+1):tempind2[2]], type='l', ylim=c(1.5,9.5), main="(d) Sunshine in Japan", xlab='Month', ylab='Sunshine (hour)', xaxt='n', cex.main=1.5, cex.lab=1.3 )
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
# 
# i <- tempcityind[1]
# lines( 1:365, Ydat2[(tempind1[i]+1):tempind1[i+1]], col=2, lwd=1.5 )
# i <- tempcityind[2]
# lines( 1:365, Ydat2[(tempind1[i]+1):tempind1[i+1]], col=3, lwd=1.5 )
# i <- tempcityind[3]
# lines( 1:365, Ydat2[(tempind1[i]+1):tempind1[i+1]], col=4, lwd=1.5 )
# legend("topright", tempcitynames, lwd=c(2,2,2), col=c(2,3,4), cex=1.4)
# 
# 
# 
# 
# i <- 1 
# dt <- data[[i]]
# 
# ## Numbering week 1. 01-01 = start of week-of-year
# dt$year <- as.factor(format(dt$date, format="%Y"))
# dt$week <- as.factor(week(dt$date))
# dt$yrwk <- as.factor(dt$year:dt$week)
# 
# # weekly data
# dt_week<-data.frame(year=as.factor(rep(unique(dt$year),each=53)),
#                     week=unique(week(dt$date)),
#                     stot=tapply(dt$stot,dt$yrwk,sum))
# 
# # correction
# tempsuic2 <- as.numeric(tapply(dt$date, as.factor(as.factor(dt$yr):as.factor(week(dt$date))), length))
# wk53ind <- which(tempsuic2 != 7)
# dt_week$stot[wk53ind] <- sapply( wk53ind, function(i){ floor(sum(dt_week$stot[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]],na.rm=T) * (7/sum(tempsuic2[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]]))) } )
# dt_week$week1 <- as.numeric(dt_week$week == 1)
# 
# # number of zero stratum
# StotStra<-tapply(dt_week$stot,dt_week$year,sum)
# drop.zero<-names(which(StotStra==0))
# 
# city_ind <- i
# NumZero[city_ind]<-length(which(dt_week$year %in% drop.zero))
# 
# # remove stratum having zero suicide count
# if(NumZero[city_ind]>0){   dt_week<-dt_week[-which(dt_week$year %in% drop.zero),]}
# 
# # Weekly data set without zero stratum
# data_week[[city_ind]]<-dt_week
# 
# plot(1:53, dt_week$stot[dt_week$year==2000], ylim=c(0,47), type='l', main="(b) Suicide counts in Japan", xlab="Month", ylab="Suicide counts", xaxt='n', cex.main=1.5, cex.lab=1.3)
# axis(1, c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5, 
#           3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), monthnames, labels = F)
# text( c(3+48/11*0, 3+48/11*1, 3+48/11*2, 3+48/11*3, 3+48/11*4, 3+48/11*5, 
#         3+48/11*6, 3+48/11*7, 3+48/11*8, 3+48/11*9, 3+48/11*10, 3+48/11*11 ), par("usr")[3] - 0.5, labels = monthnames, srt = 45, pos = 1, xpd = TRUE)
# 
# 
# collist <- c(2,3,4)
# for(i in tempcityind){
#   
#   dt <- data[[i]]
#   
#   ## Numbering week 1. 01-01 = start of week-of-year
#   dt$year <- as.factor(format(dt$date, format="%Y"))
#   dt$week <- as.factor(week(dt$date))
#   dt$yrwk <- as.factor(dt$year:dt$week)
#   
#   # weekly data
#   dt_week<-data.frame(year=as.factor(rep(unique(dt$year),each=53)),
#                       week=unique(week(dt$date)),
#                       stot=tapply(dt$stot,dt$yrwk,sum))
#   
#   # correction
#   tempsuic2 <- as.numeric(tapply(dt$date, as.factor(as.factor(dt$yr):as.factor(week(dt$date))), length))
#   wk53ind <- which(tempsuic2 != 7)
#   dt_week$stot[wk53ind] <- sapply( wk53ind, function(i){ floor(sum(dt_week$stot[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]],na.rm=T) * (7/sum(tempsuic2[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]]))) } )
#   dt_week$week1 <- as.numeric(dt_week$week == 1)
#   
#   # number of zero stratum
#   StotStra<-tapply(dt_week$stot,dt_week$year,sum)
#   drop.zero<-names(which(StotStra==0))
#   
#   city_ind <- i
#   NumZero[city_ind]<-length(which(dt_week$year %in% drop.zero))
#   
#   # remove stratum having zero suicide count
#   if(NumZero[city_ind]>0){   dt_week<-dt_week[-which(dt_week$year %in% drop.zero),]}
#   
#   # Weekly data set without zero stratum
#   data_week[[city_ind]]<-dt_week
#   
#   lines(1:53, dt_week$stot[dt_week$year==2000], type='l', col=collist[tempcityind %in% i])
# }

















# interceptvals <- NULL
# for(i in 1:n){
#   
#   dt <- data[[i]]
#   
#   ## Numbering week 1. 01-01 = start of week-of-year
#   dt$year <- as.factor(format(dt$date, format="%Y"))
#   dt$week <- as.factor(week(dt$date))
#   dt$yrwk <- as.factor(dt$year:dt$week)
#   
#   # weekly data
#   dt_week<-data.frame(year=as.factor(rep(unique(dt$year),each=53)),
#                       week=unique(week(dt$date)),
#                       stot=tapply(dt$stot,dt$yrwk,sum))
#   
#   # correction
#   tempsuic2 <- as.numeric(tapply(dt$date, as.factor(as.factor(dt$yr):as.factor(week(dt$date))), length))
#   wk53ind <- which(tempsuic2 != 7)
#   dt_week$stot[wk53ind] <- sapply( wk53ind, function(i){ floor(sum(dt_week$stot[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]],na.rm=T) * (7/sum(tempsuic2[((i-1):(i+1))[which((i-1):(i+1) %in% seq(tempsuic2))]]))) } )
#   dt_week$week1 <- as.numeric(dt_week$week == 1)
#   
#   # number of zero stratum
#   StotStra<-tapply(dt_week$stot,dt_week$year,sum)
#   drop.zero<-names(which(StotStra==0))
#   
#   city_ind <- i
#   NumZero[city_ind]<-length(which(dt_week$year %in% drop.zero))
#   
#   # remove stratum having zero suicide count
#   if(NumZero[city_ind]>0){   dt_week<-dt_week[-which(dt_week$year %in% drop.zero),]}
#   
#   # Weekly data set without zero stratum
#   data_week[[city_ind]]<-dt_week
#   
#   ### GLM
#   # seasonality basis (cyclic spline)
#   temp_spline.season <- cyclic(dt_week$week,df=4, knots=knots4)
#   
#   model <- glm(stot ~ temp_spline.season + year, data=dt_week, family=quasipoisson)
#   # COEF[[city_ind]] <- model$coefficients[2:(dfval+1)]
#   interceptvals <- c(interceptvals, coef(model)[1] + sort(coef(model)[-(1:(dfval+1))])[round((length(coef(model)) - (dfval+1))/2)])
# }
# 
# 
# mean(interceptvals[clusters == 2])
# mean(interceptvals[clusters == 3])
# mean(interceptvals[clusters == 4])
# 
# msc3[,2] <- msc3[,2] + mean(interceptvals[clusters == 2])
# msc3[,3] <- msc3[,3] + mean(interceptvals[clusters == 3])
# msc3[,4] <- msc3[,4] + mean(interceptvals[clusters == 4])
# 
# msc33 <- msc3
# plot(exp(msc33[,2]), type='l', main="(c) Expected suicide counts", ylab="Suicide counts", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5, ylim=c(5.5,11))
# axis(1, monthnum2, monthnames, labels = F)
# text(monthnum2, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ lines(exp(msc33[,clusters[which(clusters2 == i)[1]]]), col=i, lwd=3) }
# 
# 
# 
# 
# 
# par(mfrow=c(2,2))
# 
# plot(msc1[,2], type='l', ylim=c(-2,30), main="(a) Temperature", ylab="Temperature (¡ÆC)", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }
# 
# 
# plot(msc2[,2], type='l', ylim=c(2,7), main="(b) Sunshine", ylab="Sunshine (hour)", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5)
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }
# 
# msc33 <- msc3
# plot(exp(msc33[,2]), type='l', main="(c) Expected suicide counts", ylab="Suicide counts", xlab='Month', xaxt='n', cex.main=1.7, cex.lab=1.5, ylim=c(5.5,11))
# axis(1, monthnum2, monthnames, labels = F)
# text(monthnum2, par("usr")[3] - 0.2, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ lines(exp(msc33[,clusters[which(clusters2 == i)[1]]]), col=i, lwd=3) }
# 
# par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
# JapanPrefMap(clusters2, border = gray(.8), axes = TRUE)
# title(main="(d)", cex.main=1.7)
# legend(153.1, 45.56, legend = c("cluster 1", "cluster 2", "cluster 3"), col=c(2,3,4), pch = c(3,3,3), cex=1.6)
# 
# dev.off()
