################################################
source("data_generating_functions.R")
source("estimation_functions.R")

# Parameters for data generating
Mu_YX <- c(2,1,1,1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1,-1] <- ar1_cor(3, 0.3)
SigMat_YX <- SigMat_YX+0.1*diag(4)
SigMat_YX[1,1] <- 1.44
Mu_Y_T <- 1.5
Sig_Y_T <- 1.5

## True Beta
trueBetaRho <- Compute_Rho_Parameters(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)

# Sample size

n <- 500
m <- 500

# Extra Data (or CV)

n_ext <- 2000
m_ext <- 2000

data_list_extra <- Generate_Dat(n_ext, m_ext, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)

# True Distribution
yx_dist <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
coef_y_x_s_true <- c(yx_dist$Beta0, yx_dist$Beta1)
sigma_y_x_s_true <- yx_dist$VarYX

################################################

dataList <- Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
sData <- dataList$sDat
tData <- dataList$tDat
piVal <- n/(n+m)
tDat_ext <- data_list_extra$tDat
coef_y_x_s <- coef_y_x_s_true
sigma_y_x_s <- sigma_y_x_s_true
Mu_Y_S <- Mu_YX[1]
Sig_Y_S <- sqrt(SigMat_YX[1,1])

fopt <- optim(trueBetaRho, ComputeEfficientScore,
              sData = sData, tData = tData, piVal = piVal, tDat_ext = tDat_ext,
              coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, Mu_Y_S = Mu_Y_S, Sig_Y_S = Sig_Y_S)
sdVec <- ComputeCovarianceMatrix(fopt$par, sData = sData, tData = tData, piVal = piVal, tDat_ext = tDat_ext,
                                 coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, Mu_Y_S = Mu_Y_S, Sig_Y_S = Sig_Y_S)

(Cp1 <- (fopt$par[1]-1.96*sdVec[1] <= trueBetaRho[1] && fopt$par[1]+1.96*sdVec[1] >= trueBetaRho[1]))
(Cp2 <- (fopt$par[2]-1.96*sdVec[2] <= trueBetaRho[2] && fopt$par[2]+1.96*sdVec[2] >= trueBetaRho[2]))
