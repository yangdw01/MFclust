####################################################################################
# Nonparametric Bayesian latent factor model for multivariate functional data
####################################################################################
# by Daewon Yang
####################################################################################
# 2020.08.12
####################################################################################
# Canadian air pollution data
####################################################################################
# application 1 - MFLFM0
####################################################################################
rm(list=ls())

# setwd("C:\\Users\\DY\\Desktop\\MFLFM-master\\MFLFMnew\\application\\application1\\")

##########################
# load library
##########################
library(mvtnorm)
library(fda)
library(MCMCpack)
library(Rcpp)
library(RcppArmadillo)
library(maps)


##########################
# basis function
##########################
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


##########################
# data
##########################


sourceCpp("fastNEW.cpp")
load("MCCdata_Pollution_20190221.RData")
load("mcc_indicators_20190221.RData")


# mcc indicators
colnames(mcc.indicators)
colnames(mcc.indicators)[c(3,4,15,17,19,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)]

mccind_list <- c(3,4,15,17,19,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)

geo_list <- mccind_list[c(1,2)]        # all
colnames(mcc.indicators)[geo_list]

health_list <- mccind_list[c(16)]      # 08-10
colnames(mcc.indicators)[health_list]

# mccind_list[c(3,4,5,6,20)]
demo_list <- mccind_list[c(6,20)]      # 00 & 05-06, 10-11
colnames(mcc.indicators)[demo_list]

# mccind_list[c(13,14,17,18,19,21)]
socec_list <- mccind_list[c(13,14,17,19,21)]  # 10
colnames(mcc.indicators)[socec_list]

# mccind_list[c(7,8,9,10,15)]
urban_list <- mccind_list[c(7,8,9,15)] # 2000
# urban_list <- mccind_list[c(7,10)]     # 2006
colnames(mcc.indicators)[urban_list]

cat_list <- mccind_list[15]


##########################
# data preparation
##########################
funname <- c("o3","no2","pm25")
year <- 2010

cityname <-  "can8615"
citynum <- which(cities$country == cityname)
citynum_mcc <- sapply( seq(citynum), function(x){ which(mcc.indicators$city == cities[citynum[x],"city"]) } )
n <- length(citynum)

temp_list <- socec_list # socec_list  urban_list
temp_list <- temp_list[apply(mcc.indicators[citynum_mcc,temp_list], 2, function(x){ sum(is.na(x)) == 0 })]

temp_list1 <- setdiff(temp_list,cat_list)
temp_list2 <- cat_list[cat_list %in% temp_list]


# years <- year
# Ydat <- NULL
# tobs <- NULL
# nobs <- matrix(365, length(citynum), length(funname))
# for(i in citynum){
#   
#   tempdat <- dlist[[i]]
#   tempdat <- tempdat[tempdat$year %in% years,c("date","year",funname)]
#   tempdat <- tempdat[!grepl("02-29", as.character(tempdat$date)),]
#   
#   tempavg <- matrix(0, 365, length(funname))
#   temptobs <- matrix(seq(0,1,length.out = 365), 365, length(funname))
#   
#   tempavg <- tempdat[,funname]
#   
#   tempind <- !apply(tempavg,1,function(x){ sum(is.na(x))>0 })
#   
#   Ydat <- rbind(Ydat, t(t(tempavg[tempind,])))
#   tobs <- rbind(tobs, t(t(temptobs[tempind,])))
#   nobs[which(citynum %in% i),] <- sum(tempind)
# }
# 
# exind <- which(nobs < 50)
# if( length(exind)>0 ){
#   
#   temptempind <- cumsum(c(0,nobs))
#   
#   for(i in exind){
#     
#     if(nobs[i]>0){
#       
#       Ydat <- Ydat[-(temptempind[i]+1):(temptempind[i+1]),]
#       tobs <- tobs[-(temptempind[i]+1):(temptempind[i+1]),]
#       
#     }
#     
#   }
#   citynum <- citynum[-exind]
#   nobs <- t(t(nobs[-exind,]))
# }
# citynum_mcc <- sapply( seq(citynum), function(x){ which(mcc.indicators$city == cities[citynum[x],"city"]) } )
# n <- length(citynum)




years <- year
Ydat <- NULL
tobs <- NULL
exind0 <- NULL
nobs <- NULL
for(i in citynum){
  
  tempdat <- dlist[[i]]
  
  if(sum(names(tempdat) %in% funname) == length(funname)){
    
    tempdat <- tempdat[tempdat$year %in% years,c("date","year",funname)]
    tempdat <- tempdat[!grepl("02-29", as.character(tempdat$date)),]
    
    tempavg <- matrix(0, 365, length(funname))
    temptobs <- matrix(seq(0,1,length.out = 365), 365, length(funname))
    
    tempavg <- tempdat[,funname]
    
    tempind <- !apply(tempavg,1,function(x){ sum(is.na(x))>0 })
    
    Ydat <- rbind(Ydat, t(t(tempavg[tempind,])))
    tobs <- rbind(tobs, t(t(temptobs[tempind,])))
    nobs <- rbind( nobs, rep(sum(tempind), length(funname)) )
    
  }else{
    
    exind0 <- c(exind0, i)
  }
}

citynum <- citynum[!(citynum %in% exind0)]

exind <- c(which(nobs[,1] < 50))
if( length(exind)>0 ){
  
  temptempind <- cumsum(c(0,nobs[,1]))
  exindind <- NULL
  for(i in exind){ if(nobs[i,1]>0){ exindind <- c(exindind,(temptempind[i]+1):(temptempind[i+1])) } }
  citynum <- citynum[-exind]
  nobs <- t(t(nobs[-exind,]))
  
  if( length(exindind) > 0 ){
    
    Ydat <- Ydat[-exindind,]
    tobs <- tobs[-exindind,]
  }
}
citynum_mcc <- sapply( seq(citynum), function(x){ which(mcc.indicators$city == cities[citynum[x],"city"]) } )
n <- length(citynum)




Ydat[,3] <- log(Ydat[,3]+1)


############
# check
############
tempDAT <- Ydat[,3]
temptobs <- tobs[,3]
tempnobs <- nobs[,3]
tempnind <- cumsum(c(0,tempnobs))

i <- 1
plot( temptobs[(tempnind[i]+1):tempnind[i+1]], tempDAT[(tempnind[i]+1):tempnind[i+1]], type='l', ylim=quantile(tempDAT, c(0,1)) )
for(i in 1:n){ lines( temptobs[(tempnind[i]+1):tempnind[i+1]], tempDAT[(tempnind[i]+1):tempnind[i+1]] ) }

##########################
# main
##########################
niter <- 15000; 
burn <- 10000
thin <- 1
burnthin <- seq(burn+1, niter, thin)
num <- 20


basislist <- c(1,2,3,4)
Kvals <- c(65,40,20,20)

b <- 3
basis <- basislist[b]
K <- Kvals[b]

# basis
Kcs <- NULL
Bcs <- list()

for(i in 1:length(funname)){
  
  tempBasisResult <- MakeBasis(K, basis, tobs[,i], Ydat[,i])
  
  tempK <- tempBasisResult[[1]]
  tempB <- tempBasisResult[[2]]
  
  Kcs <- c(Kcs, tempK)
  Bcs[[i]] <- tempB
}


# without covariate
containCOV <- 0
containCAT <- 0
containMULT <- ifelse( length(funname)==1, 0, 1 )

Xdat <- matrix(1, 3,3)

set.seed(1300) # 1300    200 1700 2800
source("preprocess.R")

LL_result <- matrix(0, length(burnthin), n)
PS_Zmu <- list()
PS_LFLM <- list()

### MCMC iteration
for(wert in 1:niter){
  
  cat(wert, "th iteration \n")
  
  ### part 1
  tempzv = getZandV( n, M, alpha, beta, tempKpq, ppi, llamb, ZETA, a_u, b_u, a_w, b_w, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu )
  ZZ = tempzv[1:tempKpq,]
  VV = tempzv[1:tempKpq+tempKpq,]
  Eta = t( tempzv[1:n+2*tempKpq,] )
  ums = c( tempzv[2*tempKpq+n+1,] )
  mm = c( tempzv[2*tempKpq+n+2,] )
  wms = c( tempzv[2*tempKpq+n+3,] )
  Zmu = t( tempzv[1:M+2*tempKpq+n+3,] )
  
  LFLM = ZZ * VV
  r = nrow(Eta)
  
  ums = getUMS(r, ZZ, VV, a_u, b_u);
  Eta = getEta(n, r, LFLM, Sigmas, ZETA, LL, Zmu );
  
  alpha = getAlpha(ZZ, containBETA, tempKpq, beta, a_alpha, b_alpha, Hp);
  beta = getBeta(ZZ, containBETA, tempKpq, a_beta, b_beta, beta, alpha);
  
  ### part 2
  Zmu = getZmu(LL, Eta, r, wms, mm, M);
  PP = getPIS(M, eps, LL, nu);
  LL = getZis( n, M, r, PP, Eta, Zmu );
  mm = getmm(M, r, wms, Zmu);    
  wms = getwms(M, Zmu, a_w, b_w, mm, r);
  nu = getdelta( M, PP, a_nu, b_nu );
  
  ### part 3
  Sigmas = getSigmas( tempKpq, a_sigma, b_sigma, Eta, ZZ, VV, ZETA, n );
  
  if( containCAT == 1 ){ gis = getgis( tempKp, qind1, qind, qdims, Zdat, gis, n, LFLM, Eta, Sigmas, qq ); }
  if( cc>1 ){
    
    tempTHETA = getTHETA( cc, kkind, n, nobs, B, Sigmas, KKs, psi, LFLM, Eta, Ydat );
    psi = getpsi( n, a_psi, b_psi, cc, nobs, KKs, tempTHETA, Ydat, B );
    
    if(containCOV == 1){
      
      ZETA = rbind(tempTHETA, Xdat)
      
    }else{
      
      ZETA = tempTHETA
    }
    
  }else{
    
    tempTHETA = getTHETA0( n, nobs, B, Sigmas, tempK, psi0, LFLM, Eta, Ydat );
    psi0 = getpsi0( n, a_psi, b_psi, nobs, tempTHETA, Ydat, B );
    
    if(containCOV == 1){
      
      ZETA = rbind(tempTHETA, Xdat)
      
    }else{
      
      ZETA = tempTHETA
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





library(plsgenomics)

temp_ppm <- clust_result[,-1]
matrix.heatmap(temp_ppm[order(clusters),order(clusters)])



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
DIM3 <- 3
tempDAT3 <- DAT[,DIM3]
YLIM <- quantile( tempDAT, c(0,1) )
YLIM2 <- quantile( tempDAT2, c(0,1) )
YLIM3 <- quantile( tempDAT3, c(0,1) )

#### ECR algorithm (label switching)
library(abind)
library(label.switching)
tempKval <- max(LL_result)
run <- ecr(zpivot = clusters, z = LL_result, K = tempKval)

THETAmean_star <- array(NA, c(nrow(LL_result), tempK, tempKval) )
for(i in 1:nrow(LL_result)){
  
  perm <- run$permutations[i,]
  tempLL <- LL_result[i,]
  tempZmu <- PS_Zmu[[i]]
  tempZmu2 <- tempZmu[,1:tempKval]
  tempLFLM <- PS_LFLM[[i]]
  
  THETAmean_star[i,,] <- (tempLFLM %*% (tempZmu2[,perm]))[1:tempK,]
}

clusters
THETAmean_starSAMP <- THETAmean_star
THETAmean_star <- apply(THETAmean_star, c(2,3), mean)


ii <- 1
msc1 <- Bcs[[1]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_star[1:Kcs[1],]
msc2 <- Bcs[[2]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_star[1:Kcs[2]+Kcs[1],]
msc3 <- Bcs[[3]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_star[1:Kcs[3]+Kcs[1]+Kcs[2],]


msc1sample <- msc2sample <- msc3sample <- array(NA, c(nrow(LL_result), dim(msc1)[1], dim(msc1)[2]) )
msc1_lower <- msc1_upper <- msc2_lower <- msc2_upper <- msc3_lower <- msc3_upper <- msc1

for(i in 1:nrow(LL_result)){
  
  msc1sample[i,,] <- Bcs[[1]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_starSAMP[i,1:Kcs[1],]
  msc2sample[i,,] <- Bcs[[2]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_starSAMP[i,1:Kcs[2]+Kcs[1],]
  msc3sample[i,,] <- Bcs[[3]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_starSAMP[i,1:Kcs[3]+Kcs[1]+Kcs[2],]
}

for(j in 1:tempKval){
  
  msc1_lower[,j] <- c( apply(msc1sample[,,j], 2, function(x){ quantile(x,0.025, na.rm=T) }) )
  msc1_upper[,j] <- c( apply(msc1sample[,,j], 2, function(x){ quantile(x,0.975, na.rm=T) }) )
  
  msc2_lower[,j] <- c( apply(msc2sample[,,j], 2, function(x){ quantile(x,0.025, na.rm=T) }) )
  msc2_upper[,j] <- c( apply(msc2sample[,,j], 2, function(x){ quantile(x,0.975, na.rm=T) }) )
  
  msc3_lower[,j] <- c( apply(msc3sample[,,j], 2, function(x){ quantile(x,0.025, na.rm=T) }) )
  msc3_upper[,j] <- c( apply(msc3sample[,,j], 2, function(x){ quantile(x,0.975, na.rm=T) }) )
}

clusters2 <- clusters
for(w in 1:length(unique(clusters))){ clusters2[clusters == unique(clusters)[w]] <- w }



### result 
monthnum <- c(16, 47, 75, 106, 136, 167, 197,  228, 259, 289, 320, 350)
monthnames <- c("Jan.", "Feb.", "Mar.", "Apr.", "May", "June", "July", "Aug.", "Sep.", "Oct.", "Nov.", "Dec.")

# plot(msc1[,2], type='l', ylim=c(8,35), main=expression("(a) MFclust0 - O"[3]), ylab=expression("O"[3] * " (ppb)"), xlab='Month', xaxt='n')
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }
# 
# 
# plot(msc2[,2], type='l', ylim=c(5,27), main=expression("(b) MFclust0 - NO"[2]), ylab=expression("NO"[2] * " (ppb)"), xlab='Month', xaxt='n')
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }
# 
# 
# plot(msc3[,2], type='l', ylim=c(5,18), main=expression("(c) MFclust0 - PM2.5"), ylab=expression(paste("PM2.5 (", mu, "g/", m^3, ")")), xlab='Month', xaxt='n')
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 0.4, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ lines(msc3[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }
# 
# # map
# mcc_cityname <- c("Abbotsford", "Calgary", "Edmonton", "Halifax", "Hamilton", "Kingston", "Kitchener-Waterloo", "London Ontario", "Montreal", "Oakville", "Oshawa", "Ottawa", "Regina", "Sarnia",
#                   "Saint John NB", "St. John's NFL", "Sault Ste. Marie", "Saskatoon", "Thunder Bay", "Toronto", "Victoria", "Vancouver", "Windsor", "Winnipeg")
# mcc_cityname[c(2,3,16,20)]
# 
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
# map("world", "Canada" )
# points( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], pch=3, col=clusters2, lwd=2,cex=1.7 )
# # text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
# title(main="(d)", cex.main=1.5)
# legend(-53, 80, legend = c("cluster 1", "cluster 2"), col=c(1,2), pch = c(3,3,3), cex=1.5)
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[2], ((mcc.indicators[citynum_mcc,geo_list[1]]))[2]-1.2, mcc_cityname[2], cex=1, col=clusters2[2]  )
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[3], ((mcc.indicators[citynum_mcc,geo_list[1]]))[3]+2, mcc_cityname[3], cex=1, col=clusters2[3] )
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[16]-6, ((mcc.indicators[citynum_mcc,geo_list[1]]))[16]+2, mcc_cityname[16], cex=1, col=clusters2[16]  )
# segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[20], ((mcc.indicators[citynum_mcc,geo_list[1]]))[20], 
#           ((mcc.indicators[citynum_mcc,geo_list[2]]))[20]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[20]-1,col=clusters2[20]  )
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[20]+11, ((mcc.indicators[citynum_mcc,geo_list[1]]))[20]-1, mcc_cityname[20], cex=0.9, col=clusters2[20] )
























# par(mfrow=c(3,3))
# 
# plot(msc1[,2], type='l', ylim=c(8,35), main=expression("(a) MFclust0 - O"[3]), ylab=expression("O"[3] * " (ppb)"), xlab='Month', xaxt='n')
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ 
#   
#   lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
#   lines(msc1_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
#   lines(msc1_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
# }
# 
# 
# plot(msc2[,2], type='l', ylim=c(5,27), main=expression("(b) MFclust0 - NO"[2]), ylab=expression("NO"[2] * " (ppb)"), xlab='Month', xaxt='n')
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ 
#   
#   lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
#   lines(msc2_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
#   lines(msc2_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
# }
# 
# 
# plot(exp(msc3[,2])-1, type='l', ylim=c(0,17), main=expression("(c) MFclust0 - PM2.5"), ylab=expression(paste("PM2.5 (", mu, "g/", m^3, ")")), xlab='Month', xaxt='n')
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 0.4, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ 
#   
#   lines( exp(msc3[,clusters[which(clusters2 == i)[1]]])-1, col=i, lwd=3) 
#   lines( exp(msc3_upper[,clusters[which(clusters2 == i)[1]]])-1, col=i, lwd=1.5, lty=2) 
#   lines( exp(msc3_lower[,clusters[which(clusters2 == i)[1]]])-1, col=i, lwd=1.5, lty=2) 
# }
# 
# 
# plot(msc3[,2], type='l', ylim=c(0.5,3.5), main=expression("(c) MFclust0 - PM2.5"), ylab=expression(paste("PM2.5 (", mu, "g/", m^3, ")")), xlab='Month', xaxt='n')
# axis(1, monthnum, monthnames, labels = F)
# text(monthnum, par("usr")[3] - 0.4, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
# for(i in sort(unique(clusters2))){ 
#   
#   lines( msc3[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
#   lines( msc3_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
#   lines( msc3_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
# }
# 
# 
# tempDAT <- Ydat[,3]
# temptobs <- tobs[,3]
# tempnobs <- nobs[,3]
# tempnind <- cumsum(c(0,tempnobs))
# 
# i <- 1
# plot( temptobs[(tempnind[i]+1):tempnind[i+1]], tempDAT[(tempnind[i]+1):tempnind[i+1]], type='l', ylim=quantile(tempDAT, c(0,1)) )
# for(i in 1:n){ lines( temptobs[(tempnind[i]+1):tempnind[i+1]], tempDAT[(tempnind[i]+1):tempnind[i+1]], col=clusters2[i] ) }
# 
# 
# 
# tempDAT <- exp(Ydat[,3])-1
# temptobs <- tobs[,3]
# tempnobs <- nobs[,3]
# tempnind <- cumsum(c(0,tempnobs))
# 
# i <- 1
# plot( temptobs[(tempnind[i]+1):tempnind[i+1]], tempDAT[(tempnind[i]+1):tempnind[i+1]], type='l', ylim=quantile(tempDAT, c(0,1)) )
# for(i in 1:n){ lines( temptobs[(tempnind[i]+1):tempnind[i+1]], tempDAT[(tempnind[i]+1):tempnind[i+1]], col=clusters2[i] ) }
# 
# 
# 
# # map
# mcc_cityname <- c("Abbotsford", "Calgary", "Edmonton", "Halifax", "Hamilton", "Kingston", "Kitchener-Waterloo", "London Ontario", "Montreal", "Oakville", "Oshawa", "Ottawa", "Regina", "Sarnia",
#                   "Saint John NB", "St. John's NFL", "Sault Ste. Marie", "Saskatoon", "Thunder Bay", "Toronto", "Victoria", "Vancouver", "Windsor", "Winnipeg")
# mcc_cityname[c(2,3,16,20)]
# 
# par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
# map("world", "Canada" )
# points( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], pch=3, col=clusters2, lwd=2,cex=1.7 )
# # text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
# title(main="(d)", cex.main=1.5)
# legend(-53, 80, legend = c("cluster 1", "cluster 2"), col=c(1,2), pch = c(3,3,3), cex=1.5)
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[2], ((mcc.indicators[citynum_mcc,geo_list[1]]))[2]-1.2, mcc_cityname[2], cex=1, col=clusters2[2]  )
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[3], ((mcc.indicators[citynum_mcc,geo_list[1]]))[3]+2, mcc_cityname[3], cex=1, col=clusters2[3] )
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[16]-6, ((mcc.indicators[citynum_mcc,geo_list[1]]))[16]+2, mcc_cityname[16], cex=1, col=clusters2[16]  )
# segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[20], ((mcc.indicators[citynum_mcc,geo_list[1]]))[20], 
#           ((mcc.indicators[citynum_mcc,geo_list[2]]))[20]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[20]-1,col=clusters2[20]  )
# text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[20]+11, ((mcc.indicators[citynum_mcc,geo_list[1]]))[20]-1, mcc_cityname[20], cex=0.9, col=clusters2[20] )
# 
# 
# 
# 
# 
# matrix.heatmap(temp_ppm[order(clusters),order(clusters)])





plot(msc1[,2], type='l', ylim=c(8,35), main=expression("(a) MFclust0 - O"[3]), ylab=expression("O"[3] * " (ppb)"), xlab='Month', xaxt='n')
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ 
  
  lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
  lines(msc1_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
  lines(msc1_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
}


plot(msc2[,2], type='l', ylim=c(5,27), main=expression("(b) MFclust0 - NO"[2]), ylab=expression("NO"[2] * " (ppb)"), xlab='Month', xaxt='n')
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ 
  
  lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
  lines(msc2_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
  lines(msc2_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
}


plot(exp(msc3[,2])-1, type='l', ylim=c(3,17), main=expression("(c) MFclust0 - PM2.5"), ylab=expression(paste("PM2.5 (", mu, "g/", m^3, ")")), xlab='Month', xaxt='n')
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 0.4, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ 
  
  lines( exp(msc3[,clusters[which(clusters2 == i)[1]]])-1, col=i, lwd=3) 
  lines( exp(msc3_upper[,clusters[which(clusters2 == i)[1]]])-1, col=i, lwd=1.5, lty=2) 
  lines( exp(msc3_lower[,clusters[which(clusters2 == i)[1]]])-1, col=i, lwd=1.5, lty=2) 
}



# map
mcc_cityname <- c("Abbotsford", "Calgary", "Edmonton", "Halifax", "Hamilton", "Kingston", "Kitchener-Waterloo", "London Ontario", "Montreal", "Oakville", "Oshawa", "Ottawa", "Regina", "Sarnia",
                  "Saint John NB", "St. John's NFL", "Sault Ste. Marie", "Saskatoon", "Thunder Bay", "Toronto", "Victoria", "Vancouver", "Windsor", "Winnipeg")
mcc_cityname[c(2,3,16,20)]

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
map("world", "Canada" )
points( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], pch=3, col=clusters2, lwd=2,cex=1.7 )
# text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
title(main="(d)", cex.main=1.5)
legend(-63, 80, legend = c("cluster 1", "cluster 2"), col=c(1,2), pch = c(3,3,3), cex=1.5)
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[2], ((mcc.indicators[citynum_mcc,geo_list[1]]))[2]-1.2, mcc_cityname[2], cex=1, col=clusters2[2]  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[3], ((mcc.indicators[citynum_mcc,geo_list[1]]))[3]+2, mcc_cityname[3], cex=1, col=clusters2[3] )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[16]-6, ((mcc.indicators[citynum_mcc,geo_list[1]]))[16]+2, mcc_cityname[16], cex=1, col=clusters2[16]  )
segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[20], ((mcc.indicators[citynum_mcc,geo_list[1]]))[20], 
          ((mcc.indicators[citynum_mcc,geo_list[2]]))[20]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[20]-1,col=clusters2[20]  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[20]+11, ((mcc.indicators[citynum_mcc,geo_list[1]]))[20]-1, mcc_cityname[20], cex=0.9, col=clusters2[20] )





matrix.heatmap(temp_ppm[order(clusters),order(clusters)], main = "(a) MFclust0 - Heatmap" )
