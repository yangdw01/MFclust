####################################################################################
# Nonparametric Bayesian latent factor model for multivariate functional data
####################################################################################
# by Daewon Yang
####################################################################################
# 2020.04.06
####################################################################################
# simulation study 3
####################################################################################
# dynamic clustering model using dDHP
####################################################################################



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


#####################
# data
#####################
filename <- "data\\data2_ww.R"

# "data\\data2_ww.R" week     week
# "data\\data2_ws.R" week     strong
# "data\\data2_sw.R" strtong  week
# "data\\data2_ss.R" strtong  strtong

nn <- 30    # 30 50
ss <- 100  # 40 100
TT <- 4
n <- nn*TT

#####################
# source
#####################
sourceCpp("functions\\fast_dHDP_MFLFM.cpp")
source("functions\\CCR.R")


#####################
# main
#####################

set.seed(1) # 1 100 1 1

# iteration
niter <- 20000; 
burn <- 15000
thin <- 1
burnthin <- seq(burn+1, niter, thin)-1

# data generation
source(filename)

# data preprocessing
containCAT <- 0
containCOV <- 1
containBETA <- 1
containMULT <- 1

source("functions\\preprocess_dHDP.R")

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


CCR(c(clusters),truezs)


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





#######################
# estimation
#######################
cols <- c(3,2,4,6,5)

plot(tobs1[1:nobs[1,1]], msc1[,1], type='l', ylim=quantile(Ydat[,1], c(0,1)), main="(j) Case 4: Cluster-specific mean curve - dimension 1", 
     ylab=expression(y(v)), xlab=expression(v), cex.main=2, cex.lab=1.5 )
for(i in sort(unique(clusters2))){ 
  
  lines(tobs1[1:nobs[1,1]], msc1[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=3) 
  lines(tobs1[1:nobs[1,1]], msc1_upper[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
  lines(tobs1[1:nobs[1,1]], msc1_lower[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
}
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)

plot(tobs1[1:nobs[1,1]], msc2[,1], type='l', ylim=quantile(Ydat[,2], c(0,1)), main="(k) Case 4: Cluster-specific mean curve - dimension 2", 
     ylab=expression(y(v)), xlab=expression(v), cex.main=2, cex.lab=1.5)
for(i in sort(unique(clusters2))){ 
  
  lines(tobs1[1:nobs[1,1]], msc2[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=3) 
  lines(tobs1[1:nobs[1,1]], msc2_upper[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
  lines(tobs1[1:nobs[1,1]], msc2_lower[,clusters[which(clusters2 == i)[1]]], col=cols[i], lwd=1.5, lty=2) 
}
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)

matrix.heatmap(temp_ppm[order(clusters),order(clusters)], main="(l) Case 4: Heatmap", cex.main=2)

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








#######################
# covariate
#######################
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
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[4]), xlab="", main=expression("(l) Distribution of U"[4]), cex.main=1.6, cex.lab=1.3)

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
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[5]), xlab="", main=expression("(m) Distribution of U"[5]), cex.main=1.6, cex.lab=1.3)

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
boxplot(temp ~ temp2, data = temptemp, col=c(2,3,4,5,6), ylab=expression("U"[6]), xlab="", main=expression("(n) Distribution of U"[6]), cex.main=1.6, cex.lab=1.3)




#######################
# curve
#######################
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))

tt <- 1
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(a) bivariate functional data - dimension 1 at t=1", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(b) bivariate functional data - dimension 2 at t=1", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)



tt <- 2
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(c) bivariate functional data - dimension 1 at t=2", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(d) bivariate functional data - dimension 2 at t=2", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)



tt <- 3
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(e) bivariate functional data - dimension 1 at t=3", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(f) bivariate functional data - dimension 2 at t=3", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)



tt <- 4
ttind <- 1:nn + nn*(tt-1)
tt_zs <- truezs[ttind] + 1

YLIM <- quantile(Ydat11, c(0,1))
plot(te, Ydat11[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(g) bivariate functional data - dimension 1 at t=4", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat11[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)

YLIM <- quantile(Ydat22, c(0,1))
plot(te, Ydat22[,ttind[1]], type='l', ylim=YLIM, ylab=expression(y(v)), xlab=expression(v), 
     main = "(h) bivariate functional data - dimension 2 at t=4", cex.main=2.5, cex.lab=1.5)
for(i in 1:30){ lines(te, Ydat22[,ttind[i]], col=tt_zs[i]) }
legend("topleft", c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), col=c(2,3,4,5,6), lwd=c(2,2,2), cex=1.4)



