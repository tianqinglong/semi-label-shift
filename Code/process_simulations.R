library(stringr)

Extract_n_m_from_filename <- function(filename)
{
  n <- str_match(filename, "_n(.*?)_m")[2] %>% as.numeric()
  m <- str_match(filename, "_m(.*?).rds")[2] %>% as.numeric()
  
  return(list(
    n = n,
    m = m
  ))
}

Extract_Info <- function(filename, betaRhoTrue, isTrue)
{
  if (isTrue)
  {
    HatName <- "BetaRhoHat_True"
    SdName <- "BetaRhoSd_True"
  }
  else
  {
    HatName <- "BetaRhoHat_Fitted"
    SdName <- "BetaRhoSd_Fitted"
  }
  nAndm <- Extract_n_m_from_filename(filename)
  n <- nAndm$n
  m <- nAndm$m
  outList <- readRDS(paste("../SimuResults/", filename, sep = ""))
  sapply(outList, function(eachList) {
    betaRhoHat <- eachList[[HatName]]
    return(betaRhoHat)
  }) -> betaRhoHatMat
  
  betaRhoHatMat <- t(betaRhoHatMat)
  betaRhoHatMean <- colMeans(betaRhoHatMat)
  betaRhoHatBias <- betaRhoHatMean-betaRhoTrue
  betaRhoHatSE <- apply(betaRhoHatMat, MARGIN = 2, sd)
  
  sapply(outList, function(eachList) {
    betaRhoSd <- eachList[[SdName]]
    return(betaRhoSd)
  }) -> betaRhoSdMat
  betaRhoSdMat <- t(betaRhoSdMat)
  betaRhoSdMean <- colMeans(betaRhoSdMat)
  betaRho1Lower <- betaRhoHatMat[,1]-1.96*betaRhoSdMat[,1]
  betaRho1Uppoer <- betaRhoHatMat[,1]+1.96*betaRhoSdMat[,1]
  betaRho1Cp <- mean((betaRhoTrue[1] >= betaRho1Lower)*(betaRhoTrue[1] <= betaRho1Uppoer))
  
  betaRho2Lower <- betaRhoHatMat[,2]-1.96*betaRhoSdMat[,2]
  betaRho2Uppoer <- betaRhoHatMat[,2]+1.96*betaRhoSdMat[,2]
  betaRho2Cp <- mean((betaRhoTrue[2] >= betaRho2Lower)*(betaRhoTrue[2] <= betaRho2Uppoer))
  
  betaRhoCp <- c(betaRho1Cp, betaRho2Cp)
  return(
    list(
      n = n,
      m = m,
      Mean = betaRhoHatMean,
      Bias = betaRhoHatBias,
      SE = betaRhoHatSE,
      SD = betaRhoSdMean,
      CP = betaRhoCp
    )
  )
}

array_to_LaTeX <- function(arr) {
  rows <- apply(arr, MARGIN=1, paste, collapse = " & ")
  matrix_string <- paste(rows, collapse = " \\\\ ")
  return(paste("\\begin{bmatrix}", matrix_string, "\\end{bmatrix}"))
}
