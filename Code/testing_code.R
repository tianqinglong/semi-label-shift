source("data_generating_functions.R")
source("functions.R")

Mu_YX <- c(2,1,1,1)
SigMat_YX <- ar1_cor(4, 0.9)
Sig_XY[-1,-1] <- diag(3)

n <- 500
m <- 500

dat <- Generate_Y_X_Marginal(n = n, Mu_YX = Mu_YX, SigMat_YX = SigMat_YX)
XY_List <- Compute_X_Given_Y(Mu_YX = Mu_YX, SigMat_YX = SigMat_YX)

Mu_Y_T <- 1
Sig_Y_T <- 2
yTarget <- Generate_Y_Marginal_Target(m = m, Mu_Y_T = Mu_Y_T, Sig_Y_T = Sig_Y_T)
trueBetaRho <- Compute_Rho_Parameters(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)
xTarget <- Generate_X_Given_Y(yTarget, XY_List)

YX_Model_True <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
YX_Model_Fitted <- Fit_YX_Model_Source(dat)

XY_All_True <- Generate_Y_Given_X(YX_Model_True, xTarget, dat)
XY_All_Fitted <- Generate_Y_Given_X(YX_Model_Fitted, xTarget, dat)

beta_rho <- trueBetaRho
fop1 <- optim(beta_rho, ComputeEquation, yx_all = XY_All_True, sDat = dat, tDat = xTarget)
fop2 <- optim(beta_rho, ComputeEquation, yx_all = XY_All_Fitted, sDat = dat, tDat = xTarget)
trueBetaRho

ComputeCovMat(fop1$par, XY_All_True, sDat = dat, tDat = xTarget)
ComputeCovMat(fop2$par, XY_All_Fitted, sDat = dat, tDat = xTarget)
