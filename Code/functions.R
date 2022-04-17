# yx_all: is a matrix where each row corresponds to a value of x, from x_1, \dots, x_n, \dots, x_{n+m}.
# The row contains samples of y from the fitted p_s(y|x)

E_s_Rho_X <- function(beta_rho, yx_all, rho_pwr)
{
  beta_rho <- matrix(beta_rho, ncol = 1)
  rhoRowMean <- apply(yx_all, MARGIN = 1, function(y)
  {
    yTempMat <- matrix(c(y, y^2), ncol = 2, nrow = length(y))
    rhoTempVec <- exp(rho_pwr * yTempMat %*% beta_rho)
    return(mean(rhoTempVec))
  }
  )
  return(rhoRowMean)
}

# sDat: are the data from the source distribution
E_s_Rho <- function(beta_rho, sDat)
{
  beta_rho <- matrix(beta_rho, ncol = 1)
  yn <- sDat[,"Y"]
  yTempMat <- matrix(c(yn, yn^2), ncol = 2, nrow = length(yn))
  rhoVecSrc <- exp(yTempMat %*% beta_rho)
  return(mean(rhoVecSrc))
}

# c_ps: E_s(\rho)
# e_s_rho_x: E_s(\rho|\X)
# tDat: are data from the target distribution with x only
ComputeTau <- function(sDat, tDat, c_ps, e_s_rho_x, e_s_rho2_x, p1 = -1)
{
  n <- nrow(sDat)
  m <- nrow(tDat)
  if (p1 < 0)
  {
    p1 <- n/(n+m)
  }
  e_t_rho_x <- e_s_rho2_x/e_s_rho_x
  tauX <- (e_t_rho_x/p1/c_ps)/(e_t_rho_x/c_ps/p1+1/(1-p1))
  return(tauX)
}

E_t_Tau <- function(sDat, tDat, c_ps, e_s_rho_x, e_s_rho2_x, p1 = -1)
{
  n <- nrow(sDat)
  tauX <- ComputeTau(sDat, tDat, c_ps, e_s_rho_x, e_s_rho2_x, p1)
  tauXT <- tauX[-(1:n)]
  return(mean(tauXT))
}

E_t_d_log_Rho_d_Beta <- function(beta_rho, sDat, c_ps)
{
  n <- nrow(sDat)
  yn <- sDat[,"Y"]
  beta_rho <- matrix(beta_rho, ncol = 1)
  yTempMat <- matrix(c(yn, yn^2), ncol = 2, nrow = n)
  rhoVecSrc <- exp(yTempMat %*% beta_rho)
  # \beta_1: that is y
  outBeta1 <- mean(c(rhoVecSrc) * yn)/c_ps
  # \beta_2: that is y^2
  outBeta2 <- mean(c(rhoVecSrc) * (yn^2))/c_ps
  return(c(outBeta1, outBeta2))
}

E_t_d_log_Rho_d_Beta_X <- function(beta_rho, yx_all, e_s_rho_x)
{
  beta_rho <- matrix(beta_rho, ncol = 1)
  e_t_d_log_rho_d_beta_x <- matrix(0, nrow = nrow(yx_all), ncol = 2)
  for (i in 1:nrow(yx_all))
  {
    ySample <- yx_all[i,]
    yTempMat <- matrix(c(ySample, ySample^2), ncol = 2)
    rhoVecSrc <- exp(yTempMat %*% beta_rho)
    
    ySample <- matrix(ySample, ncol = 1)
    # \beta_1
    outBeta1 <- mean(rhoVecSrc*ySample)
    # \beta_2
    outBeta2 <- mean(rhoVecSrc*ySample^2)
    
    e_t_d_log_rho_d_beta_x[i,] <- c(outBeta1, outBeta2)
  }
  e_s_rho_x_mat <- matrix(e_s_rho_x, nrow = nrow(yx_all), ncol = 2)
  e_t_d_log_rho_d_beta_x <- e_t_d_log_rho_d_beta_x/e_s_rho_x_mat
  
  return(e_t_d_log_rho_d_beta_x)
}

ComputeS <- function(beta_rho, yx_all, e_s_rho_x, sDat, c_ps)
{
  e_t_d_log_rho_d_beta_x <- E_t_d_log_Rho_d_Beta_X(beta_rho, yx_all, e_s_rho_x)
  e_t_d_log_rho_d_beta <- E_t_d_log_Rho_d_Beta(beta_rho, sDat, c_ps)
  e_t_d_log_rho_d_beta <- matrix(e_t_d_log_rho_d_beta, nrow = nrow(yx_all), ncol = 2, byrow = T)
  
  return(e_t_d_log_rho_d_beta_x-e_t_d_log_rho_d_beta)
}

S_Eff_Multiplier <- function(beta_rho, n, m, yn, c_ps, p1 = -1)
{
  beta_rho <- matrix(beta_rho, ncol = 1)
  if (p1 < 0)
  {
    p1 <- n/(n+m)
  }
  yTempMat <- matrix(c(yn, yn^2), nrow = length(yn), ncol = 2)
  rhoVecSrc <- exp(yTempMat %*% beta_rho)
  scdVec1 <- c(rhoVecSrc/c_ps/p1, rep(0, m))
  scdVec2 <- c(rep(0, n), rep(-1/(1-p1), m))
  rVec <- c(rep(1, n), rep(0, m))
  tmpTerm <- rVec*scdVec1+(1-rVec)*scdVec2
  return(matrix(tmpTerm, ncol = 2, nrow = length(tmpTerm)))
}

ComputeB1 <- function(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, sDat, tDat, p1 = -1)
{
  n <- nrow(sDat)
  m <- nrow(tDat)
  if (p1 < 0)
  {
    p1 <- n/(n+m)
  }
  tauVec <- ComputeTau(sDat, tDat, c_ps, e_s_rho_x, e_s_rho2_x, p1)
  tauMat <- matrix(tauVec, ncol = 2, nrow = length(tauVec))
  sVec <- ComputeS(beta_rho, yx_all, e_s_rho_x, sDat, c_ps)
  
  e_t_tau_s <- colMeans(tauMat[-(1:n),]*sVec[-(1:n),])
  e_t_tau <- E_t_Tau(sDat, tDat, c_ps, e_s_rho_x, e_s_rho2_x, p1)
  tempTerm <- matrix(e_t_tau_s/(e_t_tau-1), ncol = 2, nrow = length(tauVec), byrow = T)
  
  return(-(1-p1)*(1-tauMat)*(sVec-tempTerm))
}

ComputeSEff <- function(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, sDat, tDat, p1 = -1)
{
  n <- nrow(sDat)
  m <- nrow(tDat)
  if (p1 < 0)
  {
    p1 <- n/(n+m)
  }
  yn <- sDat[,"Y"]
  
  s_eff_multiplier <- S_Eff_Multiplier(beta_rho, n, m, yn, c_ps, p1)
  b1 <- ComputeB1(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, sDat, tDat, p1)
  
  return(s_eff_multiplier*b1)
}

ComputeCovMatPertubation <- function(beta_rho, yx_all, sDat, tDat, ranNumExp, p1 = -1)
{
  n <- nrow(sDat)
  m <- nrow(tDat)
  if (p1 < 0)
  {
    p1 <- n/(n+m)
  }
  yn <- sDat[,"Y"]
  
  e_s_rho_x <- E_s_Rho_X(beta_rho, yx_all, 1)
  e_s_rho2_x <- E_s_Rho_X(beta_rho, yx_all, 2)
  c_ps <- E_s_Rho(beta_rho, sDat)
  
  s_eff_multiplier <- S_Eff_Multiplier(beta_rho, n, m, yn, c_ps, p1)
  b1 <- ComputeB1(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, sDat, tDat, p1)
  
  Seff <- ComputeSEff(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, sDat, tDat, p1 = -1)
  
  ranNumExp <- matrix(ranNumExp, ncol = 2, nrow = n+m)
  pertSeff <- ranNumExp*Seff
  sumPertTemp <- colMeans(pertSeff)
  
  return(sum(sumPertTemp^2))
}

ComputeEquation <- function(beta_rho, yx_all, sDat, tDat, p1 = -1)
{
  e_s_rho_x <- E_s_Rho_X(beta_rho, yx_all, 1)
  e_s_rho2_x <- E_s_Rho_X(beta_rho, yx_all, 2)
  c_ps <- E_s_Rho(beta_rho, sDat)
  sumTemp <- colMeans(ComputeSEff(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, sDat, tDat, p1))
  return(sum(sumTemp^2))
}

ComputeCovMat <- function(beta_rho, yx_all, sDat, tDat, p1 = -1)
{
  n <- nrow(sDat)
  m <- nrow(tDat)
  if (p1 < 0)
  {
    p1 <- n/(n+m)
  }
  
  e_s_rho_x <- E_s_Rho_X(beta_rho, yx_all, 1)
  e_s_rho2_x <- E_s_Rho_X(beta_rho, yx_all, 2)
  c_ps <- E_s_Rho(beta_rho, sDat)
  
  tauVec <- ComputeTau(sDat, tDat, c_ps, e_s_rho_x, e_s_rho2_x, p1)
  tauMat <- matrix(tauVec, ncol = 2, nrow = length(tauVec))
  
  sVec <- ComputeS(beta_rho, yx_all, e_s_rho_x, sDat, c_ps)
  e_t_tau_s <- colMeans(tauMat[-(1:n),]*sVec[-(1:n),])
  e_t_tau <- E_t_Tau(sDat, tDat, c_ps, e_s_rho_x, e_s_rho2_x, p1)
  tempTerm <- matrix(e_t_tau_s/(e_t_tau-1), ncol = 2, nrow = length(tauVec), byrow = T)
  tempTerm <- sVec-tempTerm
  
  tempMat <- matrix(0, nrow = length(beta_rho), ncol = length(beta_rho))
  for (i in 1:m)
  {
    tempVec <- tempTerm[(i+n),]
    tempMat <- tempMat + (1-tauVec[i+n]) * matrix(tempVec, ncol = 1) %*% matrix(tempVec, nrow = 1)
    # tempMat <- tempMat + (1-tauVec[i+n]) * matrix(tempVec, ncol = 1) %*% matrix(sVec[i+n,], nrow = 1)
  }
  tempMat <- (1-p1)*tempMat/m
  tempMat <- solve(tempMat)
  tempMat <- tempMat/(n+m)
  
  return(tempMat)
}
