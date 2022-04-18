# This function computes $E_s{\rho^{1 or 2}(y)|x}$ for the internal x, coef_y_x_s and sigma_y must be estimated using external data
E_S_RHO_X <- function(beta_rho, rho_pwr, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s, num_of_replications = 500) {
  coef_y_x_s <- matrix(coef_y_x_s, ncol = 1)
  beta_rho <- matrix(beta_rho, ncol = 1)
  
  x_mat <- cbind(1, x_mat_no_intercept)
  mean_y_x_s <- c(x_mat %*% coef_y_x_s)
  
  num_of_x <- nrow(x_mat_no_intercept)
  e_s_rho_x <- numeric(num_of_x)
  error_terms <- rnorm(num_of_replications, 0, sigma_y_x_s)
  
  for (i in 1:num_of_x) {
    y_plus_error <- mean_y_x_s[i]+error_terms
    
    yMat <- cbind(y_plus_error, y_plus_error^2)
    rho_values <- c(yMat %*% beta_rho)
    
    e_s_rho_x[i] <- mean(exp(rho_pwr*rho_values))
  }
  
  return(e_s_rho_x)
}

# This function computes $E_s(\rho(y))$, here y_s_vector is external
E_S_RHO <- function(beta_rho, y_s_external) {
  y_s_vector <- matrix(y_s_external, ncol = 1)
  y_mat <- cbind(y_s_vector, y_s_vector^2)
  
  rho_values <- c(y_mat %*% beta_rho)
  e_s_rho <- mean(exp(rho_values))
  
  return(e_s_rho)
}

# This function computes $\tau(\x)$ for the internal/external x, depending on the choice of e_s_rho_x and e_s_rho2_x
COMPUTE_TAU <- function(e_s_rho_x, e_s_rho2_x, c_ps, piVal) {
  e_t_rho_x <- e_s_rho2_x/e_s_rho_x
  tmp <- e_t_rho_x/c_ps/piVal
  tauX <- tmp/(tmp+1/(1-piVal))
  
  return(tauX)
}

# This function computes $E_t(\tau)$, x_t_mat_external.
# Here e_s_rho_x and e_s_rho2_x are related to x_t_mat_external, not the target sample
E_T_TAU <- function(x_t_mat_external, e_s_rho_x, e_s_rho2_x, c_ps, piVal) {
  tau_vec <- COMPUTE_TAU(x_t_mat_external, e_s_rho_x, e_s_rho2_x, c_ps, piVal)
  out <- mean(tau_vec)
  
  return(out)
}

# This function computes $E_t(d\log\Rho/d\beta)$, y_s_external is external
E_T_D_LOG_RHO_DIV_D_BETA <- function(beta_rho, y_s_external, c_ps) {
  beta_rho <- matrix(beta_rho, ncol = 1)
  y_s_external <- matrix(y_s_external, ncol = 1)
  yMat <- cbind(y_s_external, y_s_external^2)
  rhoVecSrc <- c(exp(yMat %*% beta_rho))
  # \beta_1: that is y
  outBeta1 <- mean(rhoVecSrc*y_s_external)/c_ps
  # \beta_2: that is y^2
  outBeta2 <- mean(rhoVecSrc*(y_s_external^2))/c_ps
  
  return(c(outBeta1, outBeta2))
}

# This function computes $\E_t(d\log\Rho/d\beta|x)$, x might be internal/external
E_T_D_LOG_RHO_DIV_D_BETA_X <- function(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s, e_s_rho_x) {
  x_mat <- cbind(1, x_mat_no_intercept)
  y_x_mean <- c(x_mat %*% coef_y_x_s)
  beta_rho <- matrix(beta_rho, ncol = 1)
  
  num_of_x <- nrow(x_mat)
  error_terms <- rnorm(num_of_x, 0, sigma_y_x_s)
  out_mat <- matrix(nrow = num_of_x, ncol = 2)
  
  for (i in 1:num_of_x) {
    y_plus_error <- y_x_mean[i]+error_terms
    y_tmp_mat <- cbind(y_plus_error, y_plus_error^2)
    rhoVecSrc <- exp(y_tmp_mat %*% beta_rho)
    
    tmp1 <- mean(rhoVecSrc*y_plus_error)
    tmp2 <- mean(rhoVecSrc*y_plus_error^2)
    
    out_mat[i,] <- c(tmp1, tmp2)
  }
  e_s_rho_x_mat <- matrix(e_s_rho_x, nrow = num_of_x, ncol = length(c(beta_rho)))
  out_mat <- out_mat/e_s_rho_x_mat
  
  return(out_mat)
}

# This function computes S(\x), for both internal/external x
Compute_S <- function(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s, e_s_rho_x, y_s_external, c_ps) {
  e_t_d_log_rho_div_d_beta_x <- E_T_D_LOG_RHO_DIV_D_BETA_X(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s, e_s_rho_x)
  e_t_d_log_rho_div_d_beta <- E_T_D_LOG_RHO_DIV_D_BETA(beta_rho, y_s_external, c_ps)
  e_t_d_log_rho_div_d_beta <- matrix(e_t_d_log_rho_div_d_beta, ncol = 2, nrow = nrow(x_mat_no_intercept), byrow = T)
  
  return(e_t_d_log_rho_div_d_beta_x-e_t_d_log_rho_div_d_beta)
}

# This function computes the efficient score function for the internal data
ComputeEfficientScore <- function(beta_rho, sData, tData, piVal, sData_ext, tDat_ext, coef_y_x_s, sigma_y_x_s) {
  num_of_source <- nrow(sData)
  num_of_target <- nrow(tData)
  num_of_total <- num_of_source+num_of_target
  
  beta_rho <- matrix(beta_rho, ncol = 1)
  y_internal <- sData[,"Y"]
  y_internal <- matrix(y_internal, ncol = 1)
  
  y_internal_mat <- cbind(y_internal, y_internal^2)
  rho_internal <- c(exp(y_internal_mat %*% beta_rho))
  
  y_s_external <- sData_ext[,"Y"]
  c_ps <- E_S_RHO(beta_rho, y_s_external)
  
  xMatAll <- rbind(sData[,-1], tData)
  
  # Find the multiplier for the S_eff
  ## For the source sample
  multiplier1 <- rho_internal/c_ps/piVal
  ## For the target sample
  multiplier2 <- rep(-1/(1-piVal), num_of_target)
  ## Concatenate
  multiplier <- matrix(c(multiplier1, multiplier2), ncol = 2, nrow = num_of_total)
  
  # Find the first part of b_1: -(1-\pi)*(1-\tau(x)) for all the internal x
  e_s_rho2_x <- E_S_RHO_X(beta_rho, 2, xMatAll, coef_y_x_s, sigma_y_x_s)
  e_s_rho_x <- E_S_RHO_X(beta_rho, 1, xMatAll, coef_y_x_s, sigma_y_x_s)
  
  tau_x_internal <- COMPUTE_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
  tau_x_internal <- matrix(c(tau_x_internal), ncol = 2, nrow = length(tau_x_internal), byrow = F)
  
  # Find the second part of b_1
  ## First find S(x), where x is internal
  s_x_internal <- Compute_S(beta_rho, xMatAll, coef_y_x_s, sigma_y_x_s, e_s_rho_x, y_s_external, c_ps)
  ## For E_t(\tau\S) and E_t(\tau) (using tDat_ext)
  e_s_rho_x_ext <- E_S_RHO_X(beta_rho, 1, tDat_ext, coef_y_x_s, sigma_y_x_s)
  e_s_rho2_x_ext <- E_S_RHO_X(beta_rho, 2, tDat_ext, coef_y_x_s, sigma_y_x_s)
  tau_x_external <- COMPUTE_TAU(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal)
  s_x_external <- Compute_S(beta_rho, tDat_ext, coef_y_x_s, sigma_y_x_s, e_s_rho_x_ext, y_s_external, c_ps)
  
  e_t_tau <- mean(tau_x_external)
  e_t_tau_s <- colMeans(matrix(tau_x_external, ncol = 2, nrow = length(tau_x_external), byrow = F)*s_x_external)
  e_t_tau_s <- matrix(e_t_tau_s, ncol = 2, nrow = num_of_total, byrow = T)
  
  b_2nd_part <- s_x_internal-e_t_tau_s/(e_t_tau-1)
  
  b1_x <- -(1-piVal)*(1-tau_x_internal)*b_2nd_part
  
  SEffMat <- multiplier*b1_x
  
  sum(colMeans(SEffMat)^2)
}

# This function computes the variance of beta hat
ComputeCovarianceMatrix <- function(beta_rho, sData, tData, piVal, sData_ext, tDat_ext, coef_y_x_s, sigma_y_x_s) {
  num_of_target <- nrow(tData)
  num_of_source <- nrow(sData)
  
  multiplier <- 1/(1-piVal)
  y_s_external <- sData_ext[,"Y"]
  c_ps <- E_S_RHO(beta_rho, y_s_external)
  
  e_s_rho2_x <- E_S_RHO_X(beta_rho, 2, tData, coef_y_x_s, sigma_y_x_s)
  e_s_rho_x <- E_S_RHO_X(beta_rho, 1, tData, coef_y_x_s, sigma_y_x_s)
  
  tau_x_t <- COMPUTE_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
  tau_x_t <- matrix(c(tau_x_t), ncol = 2, nrow = num_of_target, byrow = F)
  
  s_x_t <- Compute_S(beta_rho, tData, coef_y_x_s, sigma_y_x_s, e_s_rho_x, y_s_external, c_ps)
  
  e_s_rho_x_ext <- E_S_RHO_X(beta_rho, 1, tDat_ext, coef_y_x_s, sigma_y_x_s)
  e_s_rho2_x_ext <- E_S_RHO_X(beta_rho, 2, tDat_ext, coef_y_x_s, sigma_y_x_s)
  tau_x_external <- COMPUTE_TAU(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal)
  s_x_external <- Compute_S(beta_rho, tDat_ext, coef_y_x_s, sigma_y_x_s, e_s_rho_x_ext, y_s_external, c_ps)
  
  e_t_tau <- mean(tau_x_external)
  e_t_tau_s <- colMeans(matrix(tau_x_external, ncol = 2, nrow = length(tau_x_external), byrow = F)*s_x_external)
  e_t_tau_s <- matrix(e_t_tau_s, ncol = 2, nrow = num_of_source, byrow = T)
  
  b_2nd_part <- s_x_t-e_t_tau_s/(e_t_tau-1)
  
  outMat <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:num_of_target) {
    tmp <- matrix(b_2nd_part[i,], ncol = 1)
    outMat <- outMat+tau_x_t[i]*tmp%*%t(tmp)
  }
  outMat <- (1-piVal)*outMat/num_of_target
  cov_mat <- solve(outMat)/(num_of_target+num_of_source)
  sqrt(diag(cov_mat))
}
