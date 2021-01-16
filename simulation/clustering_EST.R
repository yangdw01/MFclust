####################################################################################
# Nonparametric Bayesian latent factor model for multivariate functional data
####################################################################################
# by Daewon Yang
####################################################################################
# 2020.04.06
####################################################################################
# simulation study 1 and 2
####################################################################################

# setwd("C:\\Users\\DY\\Desktop\\MFLFMtest\\MFLFM\\simulationNEW\\")

#####################
# load library
#####################
library(mvtnorm)
library(MCMCpack)
library(plsgenomics)
library(lpSolve)
library(corrplot)
library(fda)
library(Funclustering)
library(funHDDC)
library(maps)
library(expm)
library(R.matlab)
library(Rcpp)
library(RcppArmadillo)
library(gmfd)



#####################
# data
#####################
filename <- "data\\data1_ss.R"

# "data\\data1_ww.R" week     week
# "data\\data1_ws.R" week     strong
# "data\\data1_sw.R" strtong  week
# "data\\data1_ss.R" strtong  strtong

n <- 100    # 50 100
ss <- 100  # 50 100


#####################
# source
#####################
sourceCpp("functions\\fastNEW.cpp")
source("functions\\CCR.R")

#####################
# main
#####################

set.seed(2)
# (100,100) -> 1 4 7 6
# (100,50)  -> 1 1 5 4
# (50,100)  -> 50 8 1 39
# (50,50)   ->   34 40 2



# iteration
niter <- 15000; 
burn <- 10000
thin <- 1
burnthin <- seq(burn+1, niter, thin)-1

# data generation
source(filename)

# data preprocessing
containCAT <- 1
containCOV <- 1
containBETA <- 1
containMULT <- 1

source("functions\\preprocess.R")

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
  
  if( containCAT == 1 ){ gis = getgis( tempKp, qind1, qind, qdims, Zdatstar, gis, n, LFLM, Eta, Sigmas, qq ); }
  
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













par(mar=c(5.1, 5.1, 4.1, 2.1))

cols <- c(2,3,4,5,6)

plot(tobs1[1:nobs[1,1]], msc1[,clusters[1]], type='l', ylim=quantile(Ydat[,1], c(0,1)), main="(m) Case 4: Cluster-specific mean curve - dimension 1", 
     ylab=expression(y(v)), xlab=expression(v), cex.main=2, cex.lab=2)
for(i in sort(unique(clusters2))){ 
  
  lines(tobs1[1:nobs[1,1]], msc1[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=3) 
  lines(tobs1[1:nobs[1,1]], msc1_upper[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
  lines(tobs1[1:nobs[1,1]], msc1_lower[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
}


plot(tobs1[1:nobs[1,1]], msc2[,clusters[1]], type='l', ylim=quantile(Ydat[,2], c(0,1)), main="(n) Case 4: Cluster-specific mean curve - dimension 2", 
     ylab=expression(y(v)), xlab=expression(v), cex.main=2, cex.lab=2)
for(i in sort(unique(clusters2))){ 
  
  lines(tobs1[1:nobs[1,1]], msc2[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=3) 
  lines(tobs1[1:nobs[1,1]], msc2_upper[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
  lines(tobs1[1:nobs[1,1]], msc2_lower[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
}

plot(tobs1[1:nobs[1,1]], msc3[,clusters[1]], type='l', ylim=quantile(Ydat[,3], c(0,1)), main="(o) Case 4: Cluster-specific mean curve - dimension 3", 
     ylab=expression(y(v)), xlab=expression(v), cex.main=2, cex.lab=2)
for(i in sort(unique(clusters2))){ 
  
  lines( tobs1[1:nobs[1,1]], msc3[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=3) 
  lines( tobs1[1:nobs[1,1]], msc3_upper[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
  lines( tobs1[1:nobs[1,1]], msc3_lower[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
}

matrix.heatmap(temp_ppm[order(clusters),order(clusters)], main = "(p) Case 4: Heatmap", cex.main=2)






######################
# functional curve
######################
par(mar=c(5.1, 5.1, 4.1, 2.1))
YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(a) 3-dimensional functional data - dimension 1", cex.main=2.5, cex.lab=2 )
for(i in 1:n1){ lines(te, Ydat11[,i], col=2) }
for(i in (1:n2)+n1 ){ lines(te, Ydat11[,i],col=3) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ydat11[,i], col=4) }
for(i in (1:n4)+n1+n2+n3 ){ lines(te, Ydat11[,i], col=5) }
for(i in (1:n5)+n1+n2+n3+n4 ){ lines(te, Ydat11[,i], col=6) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2,2,2), cex=1.5)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(b) 3-dimensional functional data - dimension 2", cex.main=2.5, cex.lab=2 )
for(i in 1:n1){ lines(te, Ydat22[,i],col=2) }
for(i in (1:n2)+n1 ){ lines(te, Ydat22[,i],col=3) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ydat22[,i], col=4) }
for(i in (1:n4)+n1+n2+n3 ){ lines(te, Ydat22[,i], col=5) }
for(i in (1:n5)+n1+n2+n3+n4 ){ lines(te, Ydat22[,i], col=6) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2,2,2), cex=1.5)

YLIM <- quantile(Ydat33, c(0,1))
plot(te, Ydat33[,1], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(c) 3-dimensional functional data - dimension 3", cex.main=2.5, cex.lab=2 )
for(i in 1:n1){ lines(te, Ydat33[,i],col=2) }
for(i in (1:n2)+n1 ){ lines(te, Ydat33[,i],col=3) }
for(i in (1:n3)+n1+n2 ){ lines(te, Ydat33[,i], col=4) }
for(i in (1:n4)+n1+n2+n3 ){ lines(te, Ydat33[,i], col=5) }
for(i in (1:n5)+n1+n2+n3+n4 ){ lines(te, Ydat33[,i], col=6) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2,2,2), cex=1.5)



######################
# covariate
######################
par(mfrow=c(3,3), mar=c(5.1, 5.1, 4.1, 2.1))

i <- 1
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[1]), xlab="", 
        main=expression("(d) Distribution of U"[1]), cex.main=1.6, cex.lab=1.3)

i <- 2
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[2]), xlab="", 
        main=expression("(e) Distribution of U"[2]), cex.main=1.6, cex.lab=1.3)

i <- 3
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4), ylab=expression("U"[3]), xlab="", 
        main=expression("(f) Distribution of U"[3]), cex.main=1.6, cex.lab=1.3)

i <- 4
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[4]), xlab="", 
        main=expression("(g) Distribution of U"[4]), cex.main=1.6, cex.lab=1.3)

i <- 5
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[5]), xlab="", 
        main=expression("(h) Distribution of U"[5]), cex.main=1.6, cex.lab=1.3)

i <- 6
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[6]), xlab="", 
        main=expression("(i) Distribution of U"[6]), cex.main=1.6, cex.lab=1.3)

i <- 7
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[7]), xlab="", 
        main=expression("(j) Distribution of U"[7]), cex.main=1.6, cex.lab=1.3)

i <- 8
DIM <- i
temp <- Xdat[DIM,]
temp2 <- as.factor(c(rep("Cluster 1", n1), rep("Cluster 2", n2), rep("Cluster 3", n3), rep("Cluster 4", n4), rep("Cluster 5", n5) ))
temptemp <- data.frame(x=temp, y=temp2)
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[8]), xlab="", 
        main=expression("(k) Distribution of U"[8]), cex.main=1.6, cex.lab=1.3)


dev.off()


j <- 1
x <- t(table( Zdat[j,], rep(c(1,2,3,4,5),c(n1,n2,n3,n4,n5)) ))

colnames(x) <- c("1", "2")
barplot( x, beside=TRUE, ylab='Frequency', col=c(2,3,4,5,6), ylim=c(0,30), 
         main = expression("(l) Frequencies of z"[1]) , cex.main=2.5, cex.lab=1.6)
legend("topleft",c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), pch=c(15,15,15,15,15), col=c(2,3,4,5,6), cex=1.5 )

j <- 2
x <- t(table( Zdat[j,], rep(c(1,2,3,4,5),c(n1,n2,n3,n4,n5)) ))

colnames(x) <- c("1", "2", "3")
barplot( x, beside=TRUE, ylab='Frequency', col=c(2,3,4,5,6), ylim=c(0,30), 
         main = expression("(m) Frequencies of z"[2]) , cex.main=2.5, cex.lab=1.6)
legend("topleft",c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), pch=c(15,15,15,15,15), col=c(2,3,4,5,6), cex=1.5 )

j <- 3
x <- t(table( Zdat[j,], rep(c(1,2,3,4,5),c(n1,n2,n3,n4,n5)) ))

colnames(x) <- c("1", "2", "3")
barplot( x, beside=TRUE, ylab='Frequency', col=c(2,3,4,5,6), ylim=c(0,30), 
         main = expression("(n) Frequencies of z"[3]) , cex.main=2.5, cex.lab=1.6)
legend("topleft",c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), pch=c(15,15,15,15,15), col=c(2,3,4,5,6), cex=1.5 )





CCR(c(clusters),truezs)
