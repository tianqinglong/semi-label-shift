# Generate sample of Y given X
gen_cond_y_given_x <- function(dat, B=200)
# dat is the source distribtion sample (x,y)
{
  fout <- output_estimated_true(dat)
  sd <- fout$Sigma
  outMat <- matrix(rnorm(nrow(dat)*B, 0, sd), ncol = B)
  outMat <- outMat+matrix(fout$Fitted, ncol = B, nrow = nrow(dat))
  
  return(outMat)
}

gen_cond_y_given_t_x <- function(dat, xn, B=200)
{
  fout <- output_estimated_true(dat)
  sd <- fout$Sigma
  Coef <- fout$Coef
  fitted <- xn %*% matrix(Coef, ncol = 1)
  outMat <- matrix(rnorm(nrow(xn)*B, 0, sd), ncol = B)
  outMat <- outMat+matrix(fitted, ncol = B, nrow = nrow(xn))
  
  return(outMat)
}

# Make ingredient
E_s_rho_x <- function(betaRho, outMat)
{
  apply(outMat, 1, function(y)
  {
    yMat <- matrix(c(y, y^2), nrow = length(y)) %*% matrix(betaRho, ncol = 1)
    mean(exp(yMat))
  }  
  )
}

E_s_rho2_x <- function(betaRho, outMat)
{
  apply(outMat, 1, function(y)
  {
    yMat <- cbind(y, y^2) %*% matrix(betaRho, ncol = 1)
    mean(exp(2*yMat))
  }
  )
}

E_s_rho <- function(betaRho, dat)
# dat is the source distribtion sample (x,y)
{
  yn <- dat[,"y"]
  mean(exp(cbind(yn, yn^2) %*% matrix(betaRho, ncol = 1)))
}

Compute_tau <- function(betaRho, outMat, dat, xt, c_ps, e_s_rho_x)
# dat: source dist;
# xt: x on the target dist;
{
  n <- nrow(dat)
  m <- nrow(xt)
  p1 <- n/(n+m)
  
  e_s_rho2_x <- E_s_rho2_x(betaRho, outMat)
  
  term1 <- e_s_rho2_x/e_s_rho_x/c_ps
  tauOut <- term1/p1/(term1/p1+1/(1-p1))
  
  return(tauOut)
}

E_t_tau <- function(betaRho, outMat, dat, xt, c_ps, e_s_rho_x)
{
  n <- nrow(dat)
  tauOut <- Compute_tau(betaRho, outMat, dat, xt, c_ps, e_s_rho_x)
  return(mean(tauOut[-(1:n)]))
}

E_t_log_rho_beta <- function(betaRho, dat, c_ps)
{
  yn <- dat[,"y"]
  term0 <- exp(cbind(yn, yn^2) %*% matrix(betaRho, ncol = 1))
  term1 <- mean(term0*yn)
  term2 <- mean(term0*yn^2)
  
  return(c(term1, term2)/c_ps)
}

Compute_S <- function(betaRho, outMat, e_t_log_e_beta, e_s_rho_x)
{
  apply(outMat, 1, function(y)
  {
    mean(exp(cbind(y, y^2) %*% matrix(betaRho, ncol = 1))*y)
  }
  ) -> term1
  term1 <- term1/e_s_rho_x
  
  apply(outMat, 1, function(y)
  {
    mean(exp(cbind(y, y^2) %*% matrix(betaRho, ncol = 1))*y^2)
  }  
  ) -> term2
  term2 <- term2/e_s_rho_x
  
  out <- cbind(term1, term2)-matrix(e_t_log_e_beta, byrow = T, ncol = 2, nrow = nrow(outMat))
  
  return(out)
}


S_eff_multiplier <- function(betaRho, n, m, yn, c_ps)
{
  r <- c(rep(1, n), rep(0, m))
  p1 <- n/(n+m)
  rho_y <- exp(cbind(yn, yn^2)%*%matrix(betaRho, ncol = 1))
  out <- c(1/p1*rho_y/c_ps, rep(-1/(1-p1), m))
  out <- matrix(out, ncol = 2, nrow = n+m)
  
  return(out)
}

Generate_sampler <- function(dat, xn_t, B=200)
{
  yxs <- gen_cond_y_given_x(dat, B)
  yxt <- gen_cond_y_given_t_x(dat, xn_t, B)
  yx_all <- as.matrix(rbind(yxs, yxt))
  
  return(yx_all)
}
  
Compute_S_eff_vec <- function(betaRho, dat, xn_t, yx_all)
{
  n <- nrow(dat)
  m <- nrow(xn_t)
  p1 <- n/(n+m)
  
  c_ps <- E_s_rho(betaRho, dat)
  e_t_log_rho_beta <- E_t_log_rho_beta(betaRho, dat, c_ps)
  e_s_rho_x <- E_s_rho_x(betaRho, yx_all)
  
  tauVec <- Compute_tau(betaRho, yx_all, dat, xn_t, c_ps, e_s_rho_x)
  sVec <- Compute_S(betaRho, yx_all, e_t_log_rho_beta, e_s_rho_x)
  e_t_tau <- E_t_tau(betaRho, yx_all, dat, xn_t, c_ps, e_s_rho_x)
  
  e_t_tau_S <- c(mean(tauVec[-(1:n)]*sVec[-(1:n),1]),
                 mean(tauVec[-(1:n)]*sVec[-(1:n),2]))
  term1 <- sVec-matrix(e_t_tau_S, byrow = T, ncol = 2, nrow = n+m)/(e_t_tau-1)
  
  b1 <- -(1-p1)*matrix(1-tauVec, ncol = 2, nrow = n+m)*term1
  
  s_eff_mult <- S_eff_multiplier(betaRho, n, m, dat[,"y"], c_ps)
  
  return(s_eff_mult*b1)
}

Compute_S_eff <- function(betaRho, dat, xn_t, yx_all)
{
  sum(colMeans(Compute_S_eff_vec(betaRho, dat, xn_t, yx_all))^2)
}

Compute_beta_var <- function(betaRho, dat, xn_t, yx_all)
{
  n <- nrow(dat)
  m <- nrow(xn_t)
  p1 <- n/(n+m)
  
  c_ps <- E_s_rho(betaRho, dat)
  e_t_log_rho_beta <- E_t_log_rho_beta(betaRho, dat, c_ps)
  e_s_rho_x <- E_s_rho_x(betaRho, yx_all)
  
  tauVec <- Compute_tau(betaRho, yx_all, dat, xn_t, c_ps, e_s_rho_x)
  sVec <- Compute_S(betaRho, yx_all, e_t_log_rho_beta, e_s_rho_x)
  
  e_t_tau <- E_t_tau(betaRho, yx_all, dat, xn_t, c_ps, e_s_rho_x)
  
  e_t_tau_S <- c(
    mean(tauVec[-(1:n)]*sVec[-(1:n),1]),
    mean(tauVec[-(1:n)]*sVec[-(1:n),2])
  )
  
  term1 <- sVec-matrix(e_t_tau_S, byrow = T, ncol = 2, nrow = n+m)/(e_t_tau-1)
  SST <- matrix(0, nrow = length(betaRho), ncol = length(betaRho))
  for (i in (n+1):(n+m))
  {
    vec <- term1[i,]
    tempMat <- matrix(vec, nrow = length(betaRho), ncol = 1) %*%
      matrix(vec, nrow = 1, ncol = length(betaRho))
    tempMat <- tempMat*(1-tauVec[i])*(1-p1)
    SST <- SST+tempMat
  }
  SST <- SST/m
  
  return(SST)
}
