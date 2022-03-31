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
m <- 1000
B <- 2

MCDatList <- lapply(1:B, function(x)
{
  datList <- Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
  return(datList)
})

mclapply(MCDatList, function(datList)
{
  YX_Model_True <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
  YX_Model_Fitted <- Fit_YX_Model_Source(dat)
  
  XY_All_True <- Generate_Y_Given_X(YX_Model_True, xTarget, dat)
  XY_All_Fitted <- Generate_Y_Given_X(YX_Model_Fitted, xTarget, dat)
  
  fop1 <- optim(beta_rho, ComputeEquation, yx_all = XY_All_True, sDat = dat, tDat = xTarget)
  fop2 <- optim(beta_rho, ComputeEquation, yx_all = XY_All_Fitted, sDat = dat, tDat = xTarget)
  
  return(fop1$par)
},
mc.cores = 8
)