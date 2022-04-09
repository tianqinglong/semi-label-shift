source("data_generating_functions.R")
source("functions.R")
library(parallel)

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
B1 <- 2000
B2 <- 1000

t0 <- Sys.time()

mclapply(1:B1, function(b) {
  Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
},
mc.cores = 12) -> datList

mclapply(1:B2, function(b) {
  rexp(n+m)
},
mc.cores = 12) -> weightsList

mclapply(datList, function(datAll) {
  dat <- datAll$sDat
  xTarget <- datAll$tDat
  beta_rho <- trueBetaRho
  YX_Model_True <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
  YX_Model_Fitted <- Fit_YX_Model_Source(dat)
  
  XY_All_True <- Generate_Y_Given_X(YX_Model_True, xTarget, dat)
  XY_All_Fitted <- Generate_Y_Given_X(YX_Model_Fitted, xTarget, dat)
  
  # fop <- optim(beta_rho, ComputeEquation_cpp, yx_all = XY_All_True, sDat = dat, n = n, m = m, p1 = n/(n+m))
  # betaHat <- fop$par
  
  lapply(weightsList, function(weights) {
    fop <- optim(beta_rho, ComputeEquationPerb_cpp, yx_all = XY_All_True, sDat = dat, n = n, m = m, p1 = n/(n+m), Weights = weights)
    return(fop$par)
  }) -> beta_pertubation
  
  beta_pertubation <- matrix(unlist(beta_pertubation), byrow = T, ncol = 2)
  apply(beta_pertubation, MARGIN = 2, quantile, probs = c(0.025, 0.975)) -> CIs
  colnames(CIs) <- c("beta0", "beta1")
  
  cp1 <- (beta_rho[1] >= CIs[1,1]) & (beta_rho[1] <= CIs[2,1])
  cp2 <- (beta_rho[2] >= CIs[1,2]) & (beta_rho[2] <= CIs[2,2])
  
  return(c(cp1, cp2))
},
mc.cores = 12) -> resultsList
cp <- colMeans(matrix(unlist(resultsList), ncol = 2, byrow = T))
saveRDS(cp, file = paste("n",n,"m",m,"_cp.rds", sep = ""))

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
