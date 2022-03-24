# given yn, find max rho
max_rho <- function(yn, beta)
{
  n <- length(yn)
  ynMat <- matrix(c(yn, yn^2), ncol = 2)
  rho_list <- exp((ynMat %*% matrix(beta, ncol = 1)))
  
  return(list(cmax = max(rho_list), rho=rho_list))
}

# rejection sampling
rej_sample <- function(yn, max_out)
{
  n <- length(yn)
  
  cmax <- max_out$cmax
  rho <- max_out$rho
  rho_c <- rho/cmax

  uvec <- runif(n)
  return(uvec <= rho_c)
}

# pick selected observations
target_sample <- function(ac_list, xn, yn)
{
  return(cbind(yn[ac_list], xn[ac_list,]))
}

make_target_sample <- function(m, mu, sigma, beta, sd, beta_rho, M=5)
{
  n <- M*m
  xn <- gen_x_marginal(n, mu,sigma)
  yn <- gen_y(xn, beta, sd)
  max_out <- max_rho(yn, beta_rho)
  ac_list <- rej_sample(yn, max_out)
  
  m_temp <- sum(ac_list)
  
  while(m_temp < m)
  {
    mult <- 2*ceiling(m/m_temp)
    n <- mult*n
    
    xn <- gen_x_marginal(n, mu,sigma)
    yn <- gen_y(xn, beta, sd)
    max_out <- max_rho(yn, beta_rho)
    ac_list <- rej_sample(yn, max_out)
    
    m_temp <- sum(ac_list)
  }
  
  t_samp <- target_sample(ac_list, xn, yn)
  t_samp <- t_samp[1:m,]
  colnames(t_samp) <- c("y", "inter", "x1", "x2", "x3", "x12", "x23")
  
  return(t_samp)
}

# generate source sample
make_source_sample <- function(n, mu, sigma, beta, sd, beta_rho)
{
  xn <- gen_x_marginal(n, mu,sigma)
  yn <- gen_y(xn, beta, sd)
  
  s_samp <- cbind(yn, xn)
  colnames(s_samp) <- c("y", "inter", "x1", "x2", "x3", "x12", "x23")
  
  return(s_samp)
}
