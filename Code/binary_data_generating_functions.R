# The data generating logistic regression model
Generate_Binary_Y_Given_X <- function(xDatMatNoIntercept, f_model)
{
  n <- nrow(xDatMatNoIntercept)
  
  xDatMat <- cbind(1, xDatMatNoIntercept)
  logisCoef <- matrix(coef(f_model), ncol = 1)
  oddVec <- c(xDatMat%*%logisCoef)
  probVec <- 1/(1+exp(-oddVec))
  
  return(runif(n) <= probVec)
}

Generate_R_Given_Binary_Y <- function(yVec, gammaVec)
{
  
}

mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
mydata$rank <- factor(mydata$rank)
mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")
