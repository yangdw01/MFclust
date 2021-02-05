####################################################################################
# Nonparametric Bayesian latent factor model for multivariate functional data
####################################################################################
# by Daewon Yang
####################################################################################
# 2020.04.06
####################################################################################
# simulation study 1 and 2
####################################################################################

# setwd("C:\\Users\\DY\\Desktop\\MFLFM-master\\MFLFMnew\\simulation\\")
# setwd("C:\\Users\\DY\\Desktop\\MFLFMtest\\MFLFM\\simulationNEW\\")


Sys.time()

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
filename <- "data\\data1_ww.R"

# "data\\data1_ww.R" week     week
# "data\\data1_ws.R" week     strong
# "data\\data1_sw.R" strtong  week
# "data\\data1_ss.R" strtong  strtong

n <- 60    # 30 60
ss <- 100  # 50 100


#####################
# source
#####################
sourceCpp("functions\\fastNEW.cpp")
source("functions\\CCR.R")

#####################
# main
#####################

# iteration
niter <- 15000; 
burn <- 10000
thin <- 1
burnthin <- seq(burn+1, niter, thin)-1

NUM <- 50

MFLFM1result <- kmeans_result1 <- kmeans_result2 <- 
  MFLFM0result <- gmfd_result <- Funclust_result <- 
  kmeans_result3 <- matrix(0, NUM, n)



### functional + covariate

# case 1 - MFLFM1
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  
  # data preprocessing
  containCAT <- 1
  containCOV <- 1
  containBETA <- 1
  containMULT <- 1
  
  source("functions\\preprocess.R")
  
  # a_nu <- b_nu <- 0.5
  
  # clustering
  result <- getCLUST( burnthin, niter, eps,
                      n, M, tempKpq, tempKp, tempK, ppi, llamb, cc, containBETA, containCOV, containCAT,
                      Hp, kkind, nobs, B, KKs, Ydat, Xdat, Zdatstar, qq, qdims, qind, qind1,
                      alpha, beta, ZETA, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu, PP, nu,
                      psi, psi0, tempTHETA, gis,
                      a_u, b_u, a_w, b_w, a_alpha, b_alpha,
                      a_beta, b_beta, a_nu, b_nu, a_sigma, b_sigma, a_psi, b_psi )
  
  MFLFM1result[qwe,] <- result
}

Sys.time()


# case 3 - Kmeans
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  
  
  Ydat11 <- matrix(Ydat[[1]], length(te), n)
  Ydat22 <- matrix(Ydat[[2]], length(te), n)
  Ydat33 <- matrix(Ydat[[3]], length(te), n)
  
  CWtime <- te
  CWrange <- c(0,1)
  CWbasis <- create.bspline.basis(CWrange, nbasis=10)
  
  CWfd1 <- smooth.basisPar( CWtime, Ydat11, CWbasis, lambda=1e-2)$fd
  CWfd2 <- smooth.basisPar( CWtime, Ydat22, CWbasis, lambda=1e-2)$fd
  CWfd3 <- smooth.basisPar( CWtime, Ydat33, CWbasis, lambda=1e-2)$fd
  
  kmeans2 <- kmeans(cbind(t(rbind(CWfd1$coefs,CWfd2$coefs,CWfd3$coefs)), t(Xdat)), 3, nstart = 20)
  
  kmeans_result2[qwe,] <- kmeans2$cluster
}




### functional 

# case 1 - MFLFM0
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  
  # data preprocessing
  containCOV <- 0
  containCAT <- 0
  containBETA <- 1
  containMULT <- 1
  
  source("functions\\preprocess.R")
  
  # a_nu <- b_nu <- 0.5
  
  # clustering
  result <- getCLUST( burnthin, niter, eps,
                      n, M, tempKpq, tempKp, tempK, ppi, llamb, cc, containBETA, containCOV, containCAT,
                      Hp, kkind, nobs, B, KKs, Ydat, Xdat, Zdatstar, qq, qdims, qind, qind1,
                      alpha, beta, ZETA, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu, PP, nu,
                      psi, psi0, tempTHETA, gis,
                      a_u, b_u, a_w, b_w, a_alpha, b_alpha,
                      a_beta, b_beta, a_nu, b_nu, a_sigma, b_sigma, a_psi, b_psi )
  
  MFLFM0result[qwe,] <- result
}




# case 2 - GMFD
metrics <- c("mahalanobis","trunc","L2")

for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  tempFD <- funData( tobs[[1]][seq(nobs[[1]][1])], list(t(matrix(Ydat[[1]], ss, n)), t(matrix(Ydat[[2]], ss, n)), t(matrix(Ydat[[3]], ss, n))) )
  
  qwe
  mod1 <- gmfd_kmeans(tempFD, n.cl = 3, metric = metrics[1], p = 10^6)
  # mod3 <- gmfd_kmeans(tempFD, n.cl = 3, metric = metrics[3], p = NULL, k_trunc = NULL)
  
  gmfd_result[qwe,] <- mod1$cluster
  # gmfd_result3[qwe,] <- mod3$cluster
}

# case 3 - Funclust
for(qwe in 1:NUM){
  
  set.seed(qwe)
  source(filename)
  
  cat(qwe, "\n")
  
  # funclust & funHDDC
  Ydat1 <- matrix(Ydat[[1]], length(te), n)
  Ydat2 <- matrix(Ydat[[2]], length(te), n)
  
  # normalized
  dat <- list()
  for(i in 1:n){ dat[[i]] <- cbind( Ydat1[,i], Ydat2[,i] ) }
  
  meandat <- Reduce('+',dat)/n
  vardat <- invRdat <- list()
  for(j in 1:length(te)){
    
    temp <- matrix(0, 2, 2)
    for(i in 1:n){ temp <- temp + (dat[[i]][j,] - meandat[j,]) %*% t(dat[[i]][j,] - meandat[j,]) }
    
    vardat[[j]] <- 1/(n-1) * temp
    invRdat[[j]] <- solve( sqrtm(vardat[[j]]) )
  }
  
  for(j in 1:length(te)){
    
    invR <- invRdat[[j]]
    for(i in 1:n){ dat[[i]][j,] <- c(invR %*% dat[[i]][j,]) }
  }
  
  Ydat11 <- Ydat22 <- matrix(0, length(te), n)
  for(i in 1:n){
    
    Ydat11[,i] <- dat[[i]][,1]
    Ydat22[,i] <- dat[[i]][,2]
  }
  
  CWtime <- te
  CWrange <- c(0,1)
  CWbasis <- create.bspline.basis(CWrange, nbasis=30)
  
  CWfd1 <- smooth.basisPar( CWtime, Ydat11, CWbasis, lambda=1e-2)$fd
  CWfd2 <- smooth.basisPar( CWtime, Ydat22, CWbasis, lambda=1e-2)$fd
  CWfd <- list(CWfd1,CWfd2)
  
  res <- funclust(CWfd,K=3, thd = 0.2)
  
  Funclust_result[qwe,] <- res$cls
}



### covariate

# case 1
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  kmeans1 <- kmeans(t(Xdat), 3, nstart = 20)
  
  kmeans_result3[qwe,] <- kmeans1$cluster
}








###################
# result
###################
true_cluster <- c( rep(1,n1), rep(2, n2), rep(3,n3) )

m1_rate1 <- CCRcomplete(MFLFM1result, true_cluster)
m1_rate2 <- CCRcomplete(kmeans_result2, true_cluster)

m2_rate1 <- CCRcomplete(MFLFM0result, true_cluster)
m2_rate2 <- CCRcomplete(gmfd_result, true_cluster)
m2_rate3 <- CCRcomplete(Funclust_result, true_cluster)

m3_rate1 <- CCRcomplete(kmeans_result3, true_cluster)
m3_rate1

rate <- cbind(m1_rate1, m1_rate2, rep(NA, NUM), m2_rate1, m2_rate2, m2_rate3, rep(NA, NUM), m3_rate1)
colnames(rate) <- c("MF-LFM1", "K-means", "",  "MF-LFM0", "GMFD", "Funclust", "",  "K-means0")

boxplot(rate, main = "", ylab="Correct Classification Rate (%)", xaxt='n')
axis(side=1, at=which(!is.na(rate[1,])), labels=colnames(rate)[!is.na(rate[1,])])






Sys.time()