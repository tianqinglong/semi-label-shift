# The data generating logistic regression model
Generate_Marginal_X <- function(n, Mu_X, Sigma_X)
{
  xMat <- MASS::mvrnorm(n, Mu_X, Sigma_X)
  colnames(xMat) <- paste("X", 1:length(Mu_X), sep = "")
  
  return(xMat)
}

# Mixed: \alpha is mixed logistic coefficient
Generate_Binary_Y_Given_X <- function(xDatMatNoIntercept, alpha)
{
  n <- nrow(xDatMatNoIntercept)
  
  xDatMat <- cbind(1, xDatMatNoIntercept)
  logisCoef <- matrix(alpha, ncol = 1)
  oddVec <- c(xDatMat%*%logisCoef)
  probVec <- 1/(1+exp(-oddVec))
  
  return(as.numeric(runif(n) <= probVec))
}

Generate_R_Given_Binary_Y <- function(beta, gamma0, yVec)
{
  probVec <- 1/(1+exp(-gamma0+beta*yVec))
  n <- length(yVec)
  
  return(as.numeric(runif(n) <= probVec))
}

Generate_Source_Target_Sample <- function(numTotal, Mu_X, Sigma_X, alpha, beta, gamma0)
{
  xMarginal <- Generate_Marginal_X(numTotal, Mu_X, Sigma_X)
  yVec <- Generate_Binary_Y_Given_X(xMarginal, alpha)
  rVec <- Generate_R_Given_Binary_Y(beta, gamma0, yVec)
  
  outMat <- cbind(yVec, rVec, xMarginal)
  colnames(outMat)[1:2] <- c("Y", "R")
  
  return(outMat)
}


