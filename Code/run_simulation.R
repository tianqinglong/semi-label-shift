library(parallel)
library(Rcpp)
source("data_generating_functions.R")
source("functions.R")
sourceCpp('fast_version/functions.cpp')

Mu_YX <- c(2,1,1,1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1,-1] <- ar1_cor(3, 0.3)
SigMat_YX <- SigMat_YX+0.1*diag(4)
SigMat_YX[1,1] <- 1.44
Mu_Y_T <- 1.5
Sig_Y_T <- 1.5
trueBetaRho <- Compute_Rho_Parameters(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)

n <- 500
m <- 500

B1 <- 1800
B2 <- 1800

t0 <- Sys.time()

mclapply(1:B1, function(b) {
  Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
},
mc.cores = 12) -> datList

mclapply(1:B2, function(b) {
  rexp(n+m)
},
mc.cores = 12) -> weightsList

dataggr <- NULL
for (i in 1:length(datList)) {
  datEach <- datList[[i]]$sDat
  dataggr <- rbind(dataggr, datEach)
}

mclapply(datList, function(datAll) {
  dat <- datAll$sDat
  xTarget <- datAll$tDat
  beta_rho <- trueBetaRho
  YX_Model_True <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
  YX_Model_Fitted <- Fit_YX_Model_Source(dat)
  
  XY_All_True <- Generate_Y_Given_X(YX_Model_True, xTarget, dat)
  XY_All_Fitted <- Generate_Y_Given_X(YX_Model_Fitted, xTarget, dat)
  
  fop <- optim(beta_rho, ComputeEquation_cpp, yx_all = XY_All_True, sDat = dat, n = n, m = m, p1 = n/(n+m))
  betaHat <- fop$par
  
  covMat_cpp <- ComputeCovMat_cpp(betaHat, XY_All_True, dat, n, m, n/(n+m))
  SdVec <- sqrt(diag(covMat_cpp))
  
  c_ps_use_all <- E_s_Rho(betaHat, dataggr)
  covMat_cpp_test <- ComputeCovMat_cpp_test(betaHat, XY_All_True, dat, c_ps_use_all, n, m, n/(n+m))
  SdVec_test <- sqrt(diag(covMat_cpp_test))
  
  # lapply(weightsList, function(weights) {
  #   fop <- optim(beta_rho, ComputeEquationPerb_cpp, yx_all = XY_All_True, sDat = dat, n = n, m = m, p1 = n/(n+m), Weights = weights)
  #   return(fop$par)
  # }) -> beta_pertubation
  
  # lapply(weightsList, function(weights) {
  #   fop <- optim(beta_rho, ComputeEquationPerb_cpp, yx_all = XY_All_Fitted, sDat = dat, n = n, m = m, p1 = n/(n+m), Weights = weights)
  #   return(fop$par)
  # }) -> beta_pertubation_2

  # beta_pertubation <- matrix(unlist(beta_pertubation), byrow = T, ncol = 2)
  # apply(beta_pertubation, MARGIN = 2, quantile, probs = c(0.025, 0.975)) -> CIs
  # colnames(CIs) <- c("beta0", "beta1")
  
  # cp1 <- (beta_rho[1] >= CIs[1,1]) & (beta_rho[1] <= CIs[2,1])
  # cp2 <- (beta_rho[2] >= CIs[1,2]) & (beta_rho[2] <= CIs[2,2])
  # 
  # beta_pertubation_2 <- matrix(unlist(beta_pertubation_2), byrow = T, ncol = 2)
  # apply(beta_pertubation_2, MARGIN = 2, quantile, probs = c(0.025, 0.975)) -> CIs_2
  # colnames(CIs_2) <- c("beta0", "beta1")
  # 
  # cp3 <- (beta_rho[1] >= CIs_2[1,1]) & (beta_rho[1] <= CIs_2[2,1])
  # cp4 <- (beta_rho[2] >= CIs_2[1,2]) & (beta_rho[2] <= CIs_2[2,2])

  CI2 <- cbind(c(betaHat[1]-1.96*SdVec[1], betaHat[1]+1.96*SdVec[1]),
               c(betaHat[2]-1.96*SdVec[2], betaHat[2]+1.96*SdVec[2]))
  colnames(CI2) <- c("beta0", "beta1")
  
  CI3 <- cbind(c(betaHat[1]-1.96*SdVec_test[1], betaHat[1]+1.96*SdVec_test[1]),
               c(betaHat[2]-1.96*SdVec_test[2], betaHat[2]+1.96*SdVec_test[2]))
  colnames(CI3) <- c("beta0", "beta1")
  
  # if (CIs[1,1] <= beta_rho[1] && CIs[2,1] >= beta_rho[1]) {
  #   CP1_1 <- TRUE
  # }
  # else
  # {
  #   CP1_1 <- FALSE
  # }
  
  # if (CIs[1,2] <= beta_rho[2] && CIs[2,2] >= beta_rho[2]) {
  #   CP1_2 <- TRUE
  # }
  # else
  # {
  #   CP1_2 <- FALSE
  # }
  # CP1 <- c(CP1_1, CP1_2)
  
  if (CI2[1,1] <= beta_rho[1] && CI2[2,1] >= beta_rho[1]) {
    CP2_1 <- TRUE
  }
  else
  {
    CP2_1  <- FALSE
  }
  
  if (CI2[1,2] <= beta_rho[2] && CI2[2,2] >= beta_rho[2]) {
    CP2_2 <- TRUE
  }
  else
  {
    CP2_2  <- FALSE
  }
  CP2 <- c(CP2_1, CP2_2)
  
  if (CI3[1,1] <= beta_rho[1] && CI3[2,1] >= beta_rho[1]) {
    CP3_1 <- TRUE
  }
  else
  {
    CP3_1  <- FALSE
  }
  
  if (CI3[1,2] <= beta_rho[2] && CI3[2,2] >= beta_rho[2]) {
    CP3_2 <- TRUE
  }
  else
  {
    CP3_2  <- FALSE
  }
  CP3 <- c(CP3_1, CP3_2)
  
  return(list(True = beta_rho, Estimated = betaHat,
              # CI_1 = CIs, CP_1 = CP1,
              CI_2 = CI2, CP_2 = CP2,
              CI_3 = CI3, CP_3 = CP3))
},
mc.cores = 12) -> resultsList

# rowMeans(sapply(resultsList, function(result) {
#   result$CP_1
# }))

rowMeans(sapply(resultsList, function(result) {
  result$CP_2
}))

rowMeans(sapply(resultsList, function(result) {
  result$CP_3
}))

sapply(resultsList, function(result) {
  betaHat <- result$Estimated
  
  return(betaHat)
}) -> betaHatAll

rowMeans(betaHatAll)

sapply(resultsList, function(result) {
  len0 <- result$CI_2[2,1]-result$CI_2[1,1]
  len1 <- result$CI_2[2,2]-result$CI_2[1,2]
  
  return(c(len0, len1))
}) -> len2
rowMeans(len2)

sapply(resultsList, function(result) {
  len0 <- result$CI_3[2,1]-result$CI_3[1,1]
  len1 <- result$CI_3[2,2]-result$CI_3[1,2]
  
  return(c(len0, len1))
}) -> len3
rowMeans(len3)

# 
# cp <- colMeans(matrix(unlist(resultsList), ncol = 4, byrow = T))
# saveRDS(cp, file = paste("n",n,"m",m,"_cp.rds", sep = ""))



t1 <- Sys.time()
t1-t0

# #### Old Simulation Code
# MCDatList <- lapply(1:B, function(x)
# {
#   datList <- Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
#   return(datList)
# })
# 
# mclapply(MCDatList, function(datList)
# {
#   dat <- datList$sDat
#   xTarget <- datList$tDat
#   
#   beta_rho <- trueBetaRho
#   
#   YX_Model_True <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
#   YX_Model_Fitted <- Fit_YX_Model_Source(dat)
#   
#   XY_All_True <- Generate_Y_Given_X(YX_Model_True, xTarget, dat)
#   XY_All_Fitted <- Generate_Y_Given_X(YX_Model_Fitted, xTarget, dat)
#   
#   fop1 <- optim(beta_rho, ComputeEquation, yx_all = XY_All_True, sDat = dat, tDat = xTarget)
#   fop2 <- optim(beta_rho, ComputeEquation, yx_all = XY_All_Fitted, sDat = dat, tDat = xTarget)
#   
#   covMat_True <- ComputeCovMat(fop1$par, XY_All_True, sDat = dat, tDat = xTarget)
#   sd_True <- sqrt(diag(covMat_True))
#   
#   covMat_Fitted <- ComputeCovMat(fop2$par, XY_All_Fitted, sDat = dat, tDat = xTarget)
#   sd_Fitted <- sqrt(diag(covMat_Fitted))
#   
#   return(list(
#     BetaRhoHat_True = fop1$par,
#     BetaRhoHat_Fitted = fop2$par,
#     BetaRhoSd_True = sd_True,
#     BetaRhoSd_Fitted = sd_Fitted
#   ))
# },
# mc.cores = 12
# ) -> outList
# saveRDS(outList, file = paste("outList_n" , n, "_m", m, ".rds", sep = ""))
# 
# ##########
# 
# Extract_Info_Temp <- function(filename, betaRhoTrue, isTrue)
# {
#   if (isTrue)
#   {
#     HatName <- "BetaRhoHat_True"
#     SdName <- "BetaRhoSd_True"
#   }
#   else
#   {
#     HatName <- "BetaRhoHat_Fitted"
#     SdName <- "BetaRhoSd_Fitted"
#   }
#   nAndm <- Extract_n_m_from_filename(filename)
#   n <- nAndm$n
#   m <- nAndm$m
#   outList <- readRDS(paste("../SimuResults/", filename, sep = ""))
#   sapply(outList, function(eachList) {
#     betaRhoHat <- eachList[[HatName]]
#     return(betaRhoHat)
#   }) -> betaRhoHatMat
#   
#   betaRhoHatMat <- t(betaRhoHatMat)
#   betaRhoHatMean <- colMeans(betaRhoHatMat)
#   betaRhoHatBias <- betaRhoHatMean-betaRhoTrue
#   betaRhoHatSE <- apply(betaRhoHatMat, MARGIN = 2, sd)
#   
#   sapply(outList, function(eachList) {
#     betaRhoSd <- eachList[[SdName]]
#     return(betaRhoSd)
#   }) -> betaRhoSdMat
#   
#   betaRhoSdMat <- betaRhoSdMat/sqrt(n+m)
#   
#   betaRhoSdMat <- t(betaRhoSdMat)
#   betaRhoSdMean <- colMeans(betaRhoSdMat)
#   betaRho1Lower <- betaRhoHatMat[,1]-1.96*betaRhoSdMat[,1]
#   betaRho1Uppoer <- betaRhoHatMat[,1]+1.96*betaRhoSdMat[,1]
#   betaRho1Cp <- mean((betaRhoTrue[1] >= betaRho1Lower)*(betaRhoTrue[1] <= betaRho1Uppoer))
#   
#   betaRho2Lower <- betaRhoHatMat[,2]-1.96*betaRhoSdMat[,2]
#   betaRho2Uppoer <- betaRhoHatMat[,2]+1.96*betaRhoSdMat[,2]
#   betaRho2Cp <- mean((betaRhoTrue[2] >= betaRho2Lower)*(betaRhoTrue[2] <= betaRho2Uppoer))
#   
#   betaRhoCp <- c(betaRho1Cp, betaRho2Cp)
#   return(
#     list(
#       n = n,
#       m = m,
#       Mean = betaRhoHatMean,
#       Bias = betaRhoHatBias,
#       SE = betaRhoHatSE,
#       SD = betaRhoSdMean,
#       CP = betaRhoCp
#     )
#   )
# }
# 
# source("process_simulations.R")
# filenames <- list.files("../SimuResults/")
# 
# Extract_Info_Temp(filename = filenames[3], betaRhoTrue = trueBetaRho, isTrue = TRUE)
