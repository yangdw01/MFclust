rm(list=ls())

setwd("C:\\Users\\DY\\Desktop\\MFLFM-master\\MFLFMnew\\application\\application2\\")

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
sourceCpp("fast_dHDP_MFLFM.cpp")
load("MCCdata_Pollution_20190221.RData")
load("mcc_indicators_20190221.RData")




# mcc indicators
colnames(mcc.indicators)
colnames(mcc.indicators)[c(3,4,15,17,19,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)]

mccind_list <- c(3,4,15,17,19,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)

geo_list <- mccind_list[c(1,2)]        # all
colnames(mcc.indicators)[geo_list]


##########################
# choose country
##########################

Clist2 <- c( "bra9711",
             "can8609", "can8611",  "can8615",
             "chi9608",  "chi9615",
             "fra0010",  "fra0014",
             "ita0110",  "ita0615",
             "jap7212",
             "sa9713",
             "spa9010",  "spa9014" ,
             "tha9908",
             "uk9306Con",
             "usa7306",  "usa8506",  "usa8700"  )

ybylist_f3 <- c(7,10,16,19)
ybylist_f5 <- c(3,4,17,19)  #
ybylist_f6 <- c(2,3,4)      #
ybylist_f7 <- c(19)
ybylist_f8 <- c(19)

find <- c(5,6) # 3 5 6 7 8

##########################
# main
##########################

funname <- c("tmean","rhum","pm10","pm25","o3","no2","so2","co")[find]

cityname <-  "can8615"
citynum <- which(cities$country == cityname)
citynum_mcc <- sapply( seq(citynum), function(x){ which(mcc.indicators$city == cities[citynum[x],"city"]) } )
n <- length(citynum)

# check year
yearlist <- list()
for(i in citynum){
  
  cityi <- dlist[[i]]
  
  yearlist[[which(citynum %in% i)]] <- unique(cityi$year)
}

allyear <- yearlist[[1]]
for(i in 1:length(citynum)){ allyear <- intersect( allyear, yearlist[[i]] ) }
allyear

yearlist <- seq(allyear[1]+1, allyear[length(allyear)], 5) # seq(allyear[1], allyear[length(allyear)], 5)


all_citynum <- NULL
all_citynum_mcc <- NULL
nts <- NULL
Ydat <- NULL
tobs <- NULL
nobs <- NULL
for( years in yearlist ){
  
  year_nobs <- matrix(365, length(citynum), length(funname))
  year_Ydat <- year_tobs <- NULL
  
  for(i in citynum){
    
    tempdat <- dlist[[i]]
    tempdat <- tempdat[tempdat$year %in% years,c("date","year",funname)]
    tempdat <- tempdat[!grepl("02-29", as.character(tempdat$date)),]
    
    tempavg <- matrix(0, 365, length(funname))
    temptobs <- matrix(seq(0,1,length.out = 365), 365, length(funname))
    
    tempavg <- t(t(tempdat[,funname]))
    
    tempind <- !apply(tempavg,1,function(x){ sum(is.na(x))>0 })
    
    year_Ydat <- rbind(year_Ydat, t(t(tempavg[tempind,])))
    year_tobs <- rbind(year_tobs, t(t(temptobs[tempind,])))
    year_nobs[which(citynum %in% i),] <- sum(tempind)
  }
  
  year_exind <- which(year_nobs[,1] < 50)
  if( length(year_exind)>0 ){
    
    temptempind <- cumsum(c(0,year_nobs[,1]))
    
    for(i in year_exind){
      
      if(year_nobs[i,1]>0){
        
        year_Ydat <- year_Ydat[-(temptempind[i]+1):(temptempind[i+1]),]
        year_tobs <- year_tobs[-(temptempind[i]+1):(temptempind[i+1]),]
      }
    }
    
    year_citynum <- citynum[-year_exind]
    year_nobs <- t(t(year_nobs[-year_exind,]))
  }
  
  year_citynum_mcc <- sapply( seq(year_citynum), function(x){ which(mcc.indicators$city == cities[year_citynum[x],"city"]) } )
  year_n <- length(year_citynum)
  
  all_citynum <- c(all_citynum,year_citynum)
  all_citynum_mcc <- c(all_citynum_mcc, year_citynum_mcc)
  nts <- c(nts, year_n)
  Ydat <- rbind(Ydat, year_Ydat)
  tobs <- rbind(tobs, year_tobs)
  nobs <- rbind(nobs, year_nobs)
}
nts


### check
DAT <- Ydat

DIM <- 1
tempDAT <- DAT[,DIM]
YLIM <- quantile( tempDAT, c(0,1) )

plot( tobs[(c(0,cumsum(nobs))[1]+1):c(0,cumsum(nobs))[2]]*365, 
      tempDAT[(c(0,cumsum(nobs))[1]+1):c(0,cumsum(nobs))[2]], type='l', ylim=YLIM )
for(i in 1:sum(nts)){ 
  
  lines( tobs[(c(0,cumsum(nobs))[i]+1):c(0,cumsum(nobs))[i+1]]*365, 
         tempDAT[(c(0,cumsum(nobs))[i]+1):c(0,cumsum(nobs))[i+1]] ) 
}


### main
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

niter <- 20000; 
burn <- 10000
thin <- 1
burnthin <- seq(burn+1, niter, thin)
num <- 15

mainex <- function(){
  
  Zji_result <- matrix(0, length(burnthin), sum(nts))
  
  for(wert in 1:niter){
    
    cat(wert, "th iteration \n")
    
    ### part 1
    tempzv = getZandV( M, alpha, beta, tempKpq, ppi, llamb, ZETA, a_u, b_u, a_w, b_w, zzji, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu )
    ZZ = tempzv[1:tempKpq,]
    VV = tempzv[1:tempKpq+tempKpq,]
    Eta = t( tempzv[1:sum(nts)+2*tempKpq,] )
    ums = c( tempzv[2*tempKpq+sum(nts)+1,] )
    mm = c( tempzv[2*tempKpq+sum(nts)+2,] )
    wms = c( tempzv[2*tempKpq+sum(nts)+3,] )
    Zmu = t( tempzv[1:M+2*tempKpq+sum(nts)+3,] )
    
    LFLM = ZZ * VV
    r = nrow(Eta)
    
    ums = getUMS(r, ZZ, VV, a_u, b_u);
    Eta = getEta(r, LFLM, Sigmas, ZETA, zzji, Zmu );
    
    alpha = getAlpha(ZZ, containBETA, tempKpq, beta, a_alpha, b_alpha, Hp);
    beta = getBeta(ZZ, containBETA, tempKpq, a_beta, b_beta, beta, alpha);
    
    ### part 2
    til_wl = UPDATE__til_wl( nts, rrji, aw, bw )
    wjl = UPDATE__wjl( length(nts), til_w0, til_wl )
    til_pilk = UPDATE__til_pilk( nts, M, betak, alpha0j, rrji, zzji, eps )
    pilk = UPDATE__pilk( length(nts), M, til_pilk )
    rrji = c(UPDATE__rji( nts, wjl, pilk, zzji, t(Eta), t(Zmu) ))
    zzji = c(UPDATE__zji( nts, pilk, rrji, t(Eta), t(Zmu) ))
    alpha0j = UPDATE__alpha0j( pilk, c0, d0 )
    betak = UPDATE__betak(M, nts, eps, zzji, gam)
    gam = UPDATE__gamma( betak, gam01, gam02 )
    Zmu = getZmu(zzji, Eta, r, wms, mm, M);
    mm = getmm(M, r, wms, Zmu);
    wms = getwms(M, Zmu, a_w, b_w, mm, r);
    
    ### part 3
    Sigmas = getSigmas( tempKpq, a_sigma, b_sigma, Eta, ZZ, VV, ZETA );
    
    if( containCAT == 1 ){ gis = getgis( tempKp, qind1, qind, qdims, Zdat, gis, LFLM, Eta, Sigmas, qq ); }
    if( cc>1 ){
      
      tempTHETA = getTHETA( cc, kkind, nobs, B, Sigmas, KKs, psi, LFLM, Eta, Ydat );
      psi = getpsi( a_psi, b_psi, cc, nobs, KKs, tempTHETA, Ydat, B );
      
      if(containCOV == 1){
        
        ZETA = rbind(tempTHETA, Xdat)
        
      }else{
        
        ZETA = tempTHETA
      }
      
    }else{
      
      tempTHETA = getTHETA0( nobs, B, Sigmas, tempK, psi0, LFLM, Eta, Ydat );
      psi0 = getpsi0( a_psi, b_psi, nobs, tempTHETA, Ydat, B );
      
      if(containCOV == 1){
        
        ZETA = rbind(tempTHETA, Xdat)
        
      }else{
        
        ZETA = tempTHETA
      }
    }
    
    if( containCAT == 1 ){ ZETA = rbind( ZETA, gis ); }
    
    
    if( wert %in% burnthin ){
      
      Zji_result[which(burnthin==wert),] = zzji
    }
  }
  
  clust_result <- Dahlclust(Zji_result)
  clusters <- clust_result[,1]
  
  return(clusters)
}




timer0 <- proc.time()[3]

# # without covariate
# containCOV <- 0
# containCAT <- 0
# containMULT <- ifelse( length(funname)==1, 0, 1 )
# 
# Xdat <- matrix(1, 3,3)
# 
# set.seed(11)
# source("preprocess.R")
# tempresult1 <- mainex()

# langtidue & longitude
containCOV <- 1
containCAT <- 0
containMULT <- ifelse( length(funname)==1, 0, 1 )

cov_list <- geo_list
Xdat <- t( mcc.indicators[all_citynum_mcc,cov_list] )

set.seed(5004) # 5001 5004 5011 5019 5020 5044
source("preprocess_dHDP.R")





# tempresult2 <- mainex()
Zji_result <- matrix(0, length(burnthin), sum(nts))
PS_Zmu <- PS_LFLM <- list()

for(wert in 1:niter){
  
  cat(wert, "th iteration \n")
  
  ### part 1
  tempzv = getZandV( M, alpha, beta, tempKpq, ppi, llamb, ZETA, a_u, b_u, a_w, b_w, zzji, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu )
  ZZ = tempzv[1:tempKpq,]
  VV = tempzv[1:tempKpq+tempKpq,]
  Eta = t( tempzv[1:sum(nts)+2*tempKpq,] )
  ums = c( tempzv[2*tempKpq+sum(nts)+1,] )
  mm = c( tempzv[2*tempKpq+sum(nts)+2,] )
  wms = c( tempzv[2*tempKpq+sum(nts)+3,] )
  Zmu = t( tempzv[1:M+2*tempKpq+sum(nts)+3,] )
  
  LFLM = ZZ * VV
  r = nrow(Eta)
  
  ums = getUMS(r, ZZ, VV, a_u, b_u);
  Eta = getEta(r, LFLM, Sigmas, ZETA, zzji, Zmu );
  
  alpha = getAlpha(ZZ, containBETA, tempKpq, beta, a_alpha, b_alpha, Hp);
  beta = getBeta(ZZ, containBETA, tempKpq, a_beta, b_beta, beta, alpha);
  
  ### part 2
  til_wl = UPDATE__til_wl( nts, rrji, aw, bw )
  wjl = UPDATE__wjl( length(nts), til_w0, til_wl )
  til_pilk = UPDATE__til_pilk( nts, M, betak, alpha0j, rrji, zzji, eps )
  pilk = UPDATE__pilk( length(nts), M, til_pilk )
  rrji = c(UPDATE__rji( nts, wjl, pilk, zzji, t(Eta), t(Zmu) ))
  zzji = c(UPDATE__zji( nts, pilk, rrji, t(Eta), t(Zmu) ))
  alpha0j = UPDATE__alpha0j( pilk, c0, d0 )
  betak = UPDATE__betak(M, nts, eps, zzji, gam)
  gam = UPDATE__gamma( betak, gam01, gam02 )
  Zmu = getZmu(zzji, Eta, r, wms, mm, M);
  mm = getmm(M, r, wms, Zmu);
  wms = getwms(M, Zmu, a_w, b_w, mm, r);
  
  ### part 3
  Sigmas = getSigmas( tempKpq, a_sigma, b_sigma, Eta, ZZ, VV, ZETA );
  
  if( containCAT == 1 ){ gis = getgis( tempKp, qind1, qind, qdims, Zdat, gis, LFLM, Eta, Sigmas, qq ); }
  if( cc>1 ){
    
    tempTHETA = getTHETA( cc, kkind, nobs, B, Sigmas, KKs, psi, LFLM, Eta, Ydat );
    psi = getpsi( a_psi, b_psi, cc, nobs, KKs, tempTHETA, Ydat, B );
    
    if(containCOV == 1){
      
      ZETA = rbind(tempTHETA, Xdat)
      
    }else{
      
      ZETA = tempTHETA
    }
    
  }else{
    
    tempTHETA = getTHETA0( nobs, B, Sigmas, tempK, psi0, LFLM, Eta, Ydat );
    psi0 = getpsi0( a_psi, b_psi, nobs, tempTHETA, Ydat, B );
    
    if(containCOV == 1){
      
      ZETA = rbind(tempTHETA, Xdat)
      
    }else{
      
      ZETA = tempTHETA
    }
  }
  
  if( containCAT == 1 ){ ZETA = rbind( ZETA, gis ); }
  
  
  if( wert %in% burnthin ){
    
    Zji_result[which(burnthin==wert),] = zzji
    PS_Zmu[[which(burnthin %in% wert)]] = Zmu
    PS_LFLM[[which(burnthin %in% wert)]] = LFLM
  }
}

clust_result <- Dahlclust(Zji_result)
clusters <- clust_result[,1]
tempresult2 <- clusters





library(plsgenomics)

temp_ppm <- clust_result[,-1]
matrix.heatmap(temp_ppm[order(clusters),order(clusters)])






library(label.switching)

tempKval <- max(Zji_result)

run <- ecr( zpivot = clusters, z = Zji_result, K = tempKval )
THETAmean_star <- array(NA, c(nrow(Zji_result), tempK, tempKval) )

for(i in 1:nrow(Zji_result)){
  
  perm <- run$permutations[i,]
  tempLL <- Zji_result[i,]
  tempZmu <- PS_Zmu[[i]]
  tempZmu2 <- tempZmu[,1:tempKval]
  tempLFLM <- PS_LFLM[[i]]
  
  THETAmean_star[i,,] <- (tempLFLM %*% (tempZmu2[,perm]))[1:tempK,]
}



THETAmean_starSAMP <- THETAmean_star
THETAmean_star <- apply(THETAmean_star, c(2,3), mean)

nind <- c(0,cumsum(nobs))

ii <- 1
msc1 <- Bcs[[1]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_star[1:Kcs[1],]
msc2 <- Bcs[[2]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_star[1:Kcs[2]+Kcs[1],]



msc1sample <- msc2sample <- array(NA, c(nrow(Zji_result), dim(msc1)[1], dim(msc1)[2]) )
msc1_lower <- msc1_upper <- msc2_lower <- msc2_upper <- msc1

for(i in 1:nrow(Zji_result)){
  
  msc1sample[i,,] <- Bcs[[1]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_starSAMP[i,1:Kcs[1],]
  msc2sample[i,,] <- Bcs[[1]][(nind[ii]+1):nind[ii+1],] %*% THETAmean_starSAMP[i,1:Kcs[2]+Kcs[1],]
}

for(j in 1:tempKval){
  
  msc1_lower[,j] <- c( apply(msc1sample[,,j], 2, function(x){ quantile(x,0.025, na.rm=T) }) )
  msc1_upper[,j] <- c( apply(msc1sample[,,j], 2, function(x){ quantile(x,0.975, na.rm=T) }) )
  
  msc2_lower[,j] <- c( apply(msc2sample[,,j], 2, function(x){ quantile(x,0.025, na.rm=T) }) )
  msc2_upper[,j] <- c( apply(msc2sample[,,j], 2, function(x){ quantile(x,0.975, na.rm=T) }) )
}

clusters2 <- clusters
for(w in 1:length(unique(clusters))){ clusters2[clusters == unique(clusters)[w]] <- w }

par(mfrow=c(1,2))
plot(msc1[,1], type='l', ylim=c(0,40))
for(i in sort(unique(clusters2))){ 
  
  lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
  lines(msc1_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
  lines(msc1_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
}

plot(msc2[,1], type='l', ylim=c(0,30))
for(i in sort(unique(clusters2))){ 
  
  lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
  lines(msc2_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
  lines(msc2_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
}










clusters2 <- clusters
for(w in 1:length(unique(clusters))){ clusters2[clusters == unique(clusters)[w]] <- w }

par(mfrow=c(1,2))
plot(msc1[,1], type='l', ylim=c(0,40))
for(i in sort(unique(clusters2))){ lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i) }

plot(msc2[,1], type='l', ylim=c(0,30))
for(i in sort(unique(clusters2))){ lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i) }




monthnum <- c(16, 47, 75, 106, 136, 167, 197,  228, 259, 289, 320, 350)
monthnames <- c("Jan.", "Feb.", "Mar.", "Apr.", "May", "June", "July", "Aug.", "Sep.", "Oct.", "Nov.", "Dec.")

plot(msc1[,1], type='l', ylim=c(8,37), main=expression("(a) O"[3]*", 1987-2012"), ylab=expression("O"[3]*" (ppb)"), xlab='Month', xaxt='n')
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }
legend("topright", legend = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5", "cluster 6"), 
       col=c(1,2,4,3,5,6), lwd=c(1.5,1.5,1.5,1.5,1.5,1.5), cex=0.7)

plot(msc2[,1], type='l', ylim=c(5,25), main=expression("(b) NO"[2]*", 1987-2012"), ylab=expression("NO"[2]*" (ppb)"), xlab='Month', xaxt='n')
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) }

nts_ind <- c(0,cumsum(nts))


mcc_cityname <- c("Abbotsford", "Calgary", "Edmonton", "Halifax", "Hamilton", "Kingston", "Kitchener-Waterloo", "London Ontario", "Montreal", "Oakville", "Oshawa", "Ottawa", "Regina", "Sarnia",
                  "Saint John NB", "St. John's NFL", "Sault Ste. Marie", "Saskatoon", "Thunder Bay", "Toronto", "Victoria", "Vancouver", "Windsor", "Winnipeg")





par(mfrow=c(3,3))

i <- 1
nnind <- (nts_ind[i]+1):nts_ind[i+1]
map("world", "Canada", mar = c(4.1, 2.1, 2.1, 0.1) )
points( mcc.indicators[all_citynum_mcc[nnind],geo_list[2]], mcc.indicators[all_citynum_mcc[nnind],geo_list[1]], pch=3, col=clusters2[nnind], lwd=2, cex=2 )
# text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
title(main="(c) 1987, n=21", cex.main=1.5)
segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22], ((mcc.indicators[citynum_mcc,geo_list[1]]))[22], 
          ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1,col=2  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+10, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1.2, mcc_cityname[20], cex=1, col=2  )



i <- 2
nnind <- (nts_ind[i]+1):nts_ind[i+1]
map("world", "Canada", mar = c(4.1, 2.1, 2.1, 0.1) )
points( mcc.indicators[all_citynum_mcc[nnind],geo_list[2]], mcc.indicators[all_citynum_mcc[nnind],geo_list[1]], pch=3, col=clusters2[nnind], lwd=2, cex=2 )
# text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
title(main="(d) 1992, n=23", cex.main=1.5)
segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22], ((mcc.indicators[citynum_mcc,geo_list[1]]))[22], 
          ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1,col=2  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+10, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1.2, mcc_cityname[20], cex=1, col=2  )


i <- 3
nnind <- (nts_ind[i]+1):nts_ind[i+1]
map("world", "Canada", mar = c(4.1, 2.1, 2.1, 0.1) )
points( mcc.indicators[all_citynum_mcc[nnind],geo_list[2]], mcc.indicators[all_citynum_mcc[nnind],geo_list[1]], pch=3, col=clusters2[nnind], lwd=2, cex=2 )
# text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
title(main="(e) 1997, n=23", cex.main=1.5)
segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22], ((mcc.indicators[citynum_mcc,geo_list[1]]))[22], 
          ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1,col=2  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+10, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1.2, mcc_cityname[20], cex=1, col=2  )


i <- 4
nnind <- (nts_ind[i]+1):nts_ind[i+1]
map("world", "Canada", mar = c(4.1, 2.1, 2.1, 0.1) )
points( mcc.indicators[all_citynum_mcc[nnind],geo_list[2]], mcc.indicators[all_citynum_mcc[nnind],geo_list[1]], pch=3, col=clusters2[nnind], lwd=2, cex=2 )
# text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
title(main="(f) 2002, n=24", cex.main=1.5)
segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22], ((mcc.indicators[citynum_mcc,geo_list[1]]))[22], 
          ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1,col=5  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+10, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1.2, mcc_cityname[20], cex=1, col=5  )


i <- 5
nnind <- (nts_ind[i]+1):nts_ind[i+1]
map("world", "Canada", mar = c(4.1, 2.1, 2.1, 0.1) )
points( mcc.indicators[all_citynum_mcc[nnind],geo_list[2]], mcc.indicators[all_citynum_mcc[nnind],geo_list[1]], pch=3, col=clusters2[nnind], lwd=2, cex=2 )
# text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
title(main="(g) 2007, n=21", cex.main=1.5)
segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22], ((mcc.indicators[citynum_mcc,geo_list[1]]))[22], 
          ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1,col=6  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+10, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1.2, mcc_cityname[20], cex=1, col=6  )


i <- 6
nnind <- (nts_ind[i]+1):nts_ind[i+1]
map("world", "Canada", mar = c(4.1, 2.1, 2.1, 0.1) )
points( mcc.indicators[all_citynum_mcc[nnind],geo_list[2]], mcc.indicators[all_citynum_mcc[nnind],geo_list[1]], pch=3, col=clusters2[nnind], lwd=2, cex=2 )
# text( mcc.indicators[citynum_mcc,geo_list[2]], mcc.indicators[citynum_mcc,geo_list[1]], 1:n, cex=0.8, col=tempresult )
title(main="(h) 2012, n=24", cex.main=1.5)
segments( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22], ((mcc.indicators[citynum_mcc,geo_list[1]]))[22], 
          ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+5, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1,col=6  )
text( ((mcc.indicators[citynum_mcc,geo_list[2]]))[22]+10, ((mcc.indicators[citynum_mcc,geo_list[1]]))[22]-1.2, mcc_cityname[20], cex=1, col=6  )



plot(msc1[,1], type='l', ylim=c(8,37), main=expression("(a) O"[3]*", 1987-2012"), ylab=expression("O"[3]*" (ppb)"), xlab='Month', xaxt='n')
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ 
  
  lines(msc1[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
  lines(msc1_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
  lines(msc1_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
}
legend("topright", legend = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5", "cluster 6"), 
       col=c(1,2,4,3,5,6), lwd=c(1.5,1.5,1.5,1.5,1.5,1.5), cex=0.7)

plot(msc2[,1], type='l', ylim=c(5,25), main=expression("(b) NO"[2]*", 1987-2012"), ylab=expression("NO"[2]*" (ppb)"), xlab='Month', xaxt='n')
axis(1, monthnum, monthnames, labels = F)
text(monthnum, par("usr")[3] - 1, labels = monthnames, srt = 45, pos = 1, xpd = TRUE, cex=0.85)
for(i in sort(unique(clusters2))){ 
  
  lines(msc2[,clusters[which(clusters2 == i)[1]]], col=i, lwd=3) 
  lines(msc2_upper[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
  lines(msc2_lower[,clusters[which(clusters2 == i)[1]]], col=i, lwd=1.5, lty=2) 
}

matrix.heatmap(temp_ppm[order(clusters),order(clusters)])

length( unique(clusters) )
