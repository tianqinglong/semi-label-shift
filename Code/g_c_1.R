# Generate data from setting 1

gen_x_marginal <- function(n, mu, sigma)
{
  xn <- MASS::mvrnorm(n, mu, sigma)

  x12 <- xn[,1]*xn[,2]
  x23 <- xn[,2]*xn[,3]
  
  xn <- cbind(1, xn, x12, x23)
  return(xn)
}

gen_y <- function(xn, beta, sd)
{
  n <- nrow(xn)
  y <- xn %*% as.matrix(beta, ncol = 1) + as.matrix(rnorm(n, 0, sd), ncol = 1)
  
  return(y)
}
