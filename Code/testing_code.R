source("data_generating_functions.R")
source("functions.R")
Rcpp::sourceCpp("fast_version/functions.cpp")

Mu_YX <- c(2,1,1,1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1,-1] <- ar1_cor(3, 0.5)
SigMat_YX <- SigMat_YX+diag(4)*0.2

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
fop1
fop1_cpp <- optim(beta_rho, ComputeEquation_cpp, yx_all = XY_All_True, sDat = dat, n = n, m = m, p1 = n/(n+m))
fop1_cpp

covMat_cpp <- ComputeCovMat_cpp(beta_rho, XY_All_True, dat, n, m, n/(n+m))
covMat <- ComputeCovMat(beta_rho = beta_rho, yx_all = XY_All_True, sDat = dat, tDat = xTarget, p1 = n/(n+m))

sqrt(diag(covMat_cpp))

##

library(parallel)

B <- 1000
weightsList <- mclapply(1:B, function(b) {
  return(rexp(n+m))
},
mc.cores = 12
)

mclapply(weightsList, function(weights) {
  fop <- optim(beta_rho, ComputeEquationPerb_cpp, yx_all = XY_All_True, sDat = dat, n = n, m = m, p1 = n/(n+m), Weights = weights)
  return(fop$par)
},
mc.cores = 12
) -> beta_pertubation

beta_pertubation <- matrix(unlist(beta_pertubation), byrow = T, ncol = 2)
apply(beta_pertubation, MARGIN = 2, quantile, probs = c(0.025, 0.975))

##

nm <- 300
Mu_X <- rep(1,3)
Sigma_X <- diag(3)
Alpha <- c(-1, 2, -2, 2)
beta <- -1
gamma0 <- 0

datTarSrc <- Generate_Source_Target_Sample(nm, Mu_X, Sigma_X, Alpha, beta, gamma0)
datTarSrc <- as.data.frame(datTarSrc)
table(Y = datTarSrc[,"Y"], R = datTarSrc[,"R"])

