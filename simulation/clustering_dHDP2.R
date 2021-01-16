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
# For DP models
####################################################################################

# setwd("C:\\Users\\DY\\Desktop\\MFLFM-master\\MFLFMnew\\simulation\\")

timer0 <- proc.time()[3]

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
filename <- "data\\data2_ww.R"

# "data\\data2_ww.R" week     week
# "data\\data2_ws.R" week     strong
# "data\\data2_sw.R" strtong  week
# "data\\data2_ss.R" strtong  strtong

nn <- 30    # 30 50
ss <- 100  # 40 100
TT <- 4

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

NUM <- 100 
n <- nn*TT

dHDP_MFLFM1result <- MFLFM1result <- kmeans_result1 <- kmeans_result2 <- 
  dHDP_MFLFM0result <- MFLFM0result <- gmfd_result <- Funclust_result <- 
  kmeans_result3 <- matrix(0, NUM, n)



# 6:21


### functional + covariate

# # case 1 - dHDP-MFLFM1
# for(qwe in 1:NUM){
# 
#   set.seed(qwe)
# 
#   # data generation
#   source(filename)
# 
#   # data preprocessing
#   containCOV <- 1
#   containCAT <- 0
#   containBETA <- 1
#   containMULT <- 1
# 
#   source("functions\\preprocess_dHDP.R")
# 
#   # clustering
#   result <- mainex()
# 
#   # result
#   dHDP_MFLFM1result[qwe,] <- result
# }


# case 2 - MFLFM1
for(qwe in 1:NUM){

  set.seed(qwe)

  # data generation
  source(filename)

  # data preprocessing
  containCOV <- 1
  containCAT <- 0
  containBETA <- 1
  containMULT <- 1

  source("functions\\preprocess.R")

  # clustering
  result <- getCLUST( burnthin, niter, eps,
                      n, M, tempKpq, tempKp, tempK, ppi, llamb, cc, containBETA, containCOV, containCAT,
                      Hp, kkind, nobs, B, KKs, Ydat, Xdat, Zdatstar, qq, qdims, qind, qind1,
                      alpha, beta, ZETA, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu, PP, nu,
                      psi, psi0, tempTHETA, gis,
                      a_u, b_u, a_w, b_w, a_alpha, b_alpha,
                      a_beta, b_beta, a_nu, b_nu, a_sigma, b_sigma, a_psi, b_psi )

  # result
  MFLFM1result[qwe,] <- result
}

# case 3 - Kmeans
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  
  tempcoef <- matrix( NA, n, sum(sapply(Bcs, function(x){ dim(x)[2] })) )
  tempnobs <- nobs[[1]]
  tempnobs_ind <- cumsum(c(0,tempnobs))
  
  for(i in 1:n){
    
    temptempcoef <- NULL
    tempiind <- (tempnobs_ind[i]+1):tempnobs_ind[i+1]
    for(ii in 1:cc){ temptempcoef <- c(temptempcoef, lm( Ydat[[ii]][tempiind] ~ 0 + Bcs[[ii]][tempiind,] )$coefficients) }
    
    tempcoef[i,] <- temptempcoef
  }
  
  kmeans2 <- kmeans(cbind(tempcoef, t(Xdat)), 5, nstart = 20)
  
  kmeans_result2[qwe,] <- kmeans2$cluster
}

### functional

# # case 1 - dHDP-MFLFM0
# for(qwe in 1:NUM){
#   
#   set.seed(qwe)
#   
#   # data generation
#   source(filename)
#   
#   # data preprocessing
#   containCOV <- 0
#   containCAT <- 0
#   containBETA <- 1
#   containMULT <- 1
#   
#   source("functions\\preprocess_dHDP.R")
#   
#   # clustering
#   result <- mainex()
#   
#   # result
#   dHDP_MFLFM0result[qwe,] <- result
# }


# case 2 - MFLFM0
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
  
  # clustering
  # clustering
  result <- getCLUST( burnthin, niter, eps,
                      n, M, tempKpq, tempKp, tempK, ppi, llamb, cc, containBETA, containCOV, containCAT,
                      Hp, kkind, nobs, B, KKs, Ydat, Xdat, Zdatstar, qq, qdims, qind, qind1,
                      alpha, beta, ZETA, LL, ZZ, VV, Eta, ums, Sigmas, mm, wms, Zmu, PP, nu,
                      psi, psi0, tempTHETA, gis,
                      a_u, b_u, a_w, b_w, a_alpha, b_alpha,
                      a_beta, b_beta, a_nu, b_nu, a_sigma, b_sigma, a_psi, b_psi )
  
  # result
  MFLFM0result[qwe,] <- result
}

contime <- (proc.time()[3] - timer0)/60

# case 3 - GMFD
metrics <- c("mahalanobis","trunc","L2")

for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  tempFD <- funData( tobs[[1]][seq(nobs[[1]][1])], list(t(matrix(Ydat[[1]], ss, n)), t(matrix(Ydat[[2]], ss, n))) )
  
  try({
    
    mod1 <- gmfd_kmeans(tempFD, n.cl = 5, metric = metrics[1], p = 10^6)
    gmfd_result[qwe,] <- mod1$cluster
    
  }, silent=TRUE)
}

# case 4 - Funclust
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
  
  res <- funclust(CWfd,K=5, thd = 0.2)
  
  Funclust_result[qwe,] <- res$cls
}


### covariate

# case 1 - K means
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  kmeans1 <- kmeans(t(Xdat), 5, nstart = 10)
  
  kmeans_result3[qwe,] <- kmeans1$cluster
}










# ###################
# # result
# ###################
# true_clusters <- matrix(NA, NUM, n)
# for(qwe in 1:NUM){
#   
#   set.seed(qwe)
#   
#   # data generation
#   source(filename)
#   true_clusters[qwe,] <- truezs
# }
# 
# m1_rate1 <- sapply(1:NUM, function(x){ CCR( dHDP_MFLFM1result[x,], true_clusters[x,] ) })
# m1_rate2 <- sapply(1:NUM, function(x){ CCR( MFLFM1result[x,], true_clusters[x,] ) })
# m1_rate3 <- sapply(1:NUM, function(x){ CCR( kmeans_result2[x,], true_clusters[x,] ) })
# 
# m2_rate1 <- sapply(1:NUM, function(x){ CCR( dHDP_MFLFM0result[x,], true_clusters[x,] ) })
# m2_rate2 <- sapply(1:NUM, function(x){ CCR( MFLFM0result[x,], true_clusters[x,] ) })
# m2_rate3 <- sapply(1:NUM, function(x){ CCR( gmfd_result[x,], true_clusters[x,] ) })
# m2_rate4 <- sapply(1:NUM, function(x){ CCR( Funclust_result[x,], true_clusters[x,] ) })
# 
# m3_rate1 <- sapply(1:NUM, function(x){ CCR( kmeans_result3[x,], true_clusters[x,] ) })
# 
# rate <- cbind(m1_rate1, m1_rate2, m1_rate3, rep(NA, NUM), 
#               m2_rate1, m2_rate2, m2_rate3, m2_rate4, rep(NA, NUM), m3_rate1)
# colnames(rate) <- c("dHDP-MF-LFM1", "MF-LFM1", "K-means", "",  
#                     "dHDP-MF-LFM0", "MF-LFM0", "GMFD", "Funclust", "",  
#                     "K-means0")
# 
# boxplot(rate, main = "", ylab="Correct Classification Rate (%)", xaxt='n')
# axis(side=1, at=which(!is.na(rate[1,])), labels=colnames(rate)[!is.na(rate[1,])])




###################
# result
###################
true_clusters <- matrix(NA, NUM, n)
for(qwe in 1:NUM){
  
  set.seed(qwe)
  
  # data generation
  source(filename)
  true_clusters[qwe,] <- truezs
}

# m1_rate1 <- sapply(1:NUM, function(x){ CCR( dHDP_MFLFM1result[x,], true_clusters[x,] ) })
m1_rate2 <- sapply(1:NUM, function(x){ CCR( MFLFM1result[x,], true_clusters[x,] ) })
m1_rate3 <- sapply(1:NUM, function(x){ CCR( kmeans_result2[x,], true_clusters[x,] ) })

# m2_rate1 <- sapply(1:NUM, function(x){ CCR( dHDP_MFLFM0result[x,], true_clusters[x,] ) })
m2_rate2 <- sapply(1:NUM, function(x){ CCR( MFLFM0result[x,], true_clusters[x,] ) })
m2_rate3 <- sapply(1:NUM, function(x){ CCR( gmfd_result[x,], true_clusters[x,] ) })
m2_rate4 <- sapply(1:NUM, function(x){ CCR( Funclust_result[x,], true_clusters[x,] ) })

m3_rate1 <- sapply(1:NUM, function(x){ CCR( kmeans_result3[x,], true_clusters[x,] ) })

rate <- cbind(m1_rate2, m1_rate3, rep(NA, NUM), 
              m2_rate2, m2_rate3, m2_rate4, rep(NA, NUM), m3_rate1)
colnames(rate) <- c("MF-LFM1", "K-means", "",  
                    "MF-LFM0", "GMFD", "Funclust", "",  
                    "K-means0")

boxplot(rate, main = "", ylab="Correct Classification Rate (%)", xaxt='n')
axis(side=1, at=which(!is.na(rate[1,])), labels=colnames(rate)[!is.na(rate[1,])])





# save(MFLFM1result, file="MFLFM1result.RData")
# save(MFLFM0result, file="MFLFM0result.RData")
# save(kmeans_result2, file="kmeans_result2.RData")
# save(kmeans_result3, file="kmeans_result3.RData")
# save(Funclust_result, file="Funclust_result.RData")
# save(gmfd_result, file="gmfd_result.RData")