# Fit the correct model
output_estimated_true <- function(dat)
{
  f1 <- lm(y~.-inter, data = as.data.frame(dat))
  
  return(list(Fitted=f1$fitted.values, Sigma=sigma(f1), Coef = f1$coefficients))
}

output_estimated_wrong <- function(dat)
{
  f1 <- lm(y~x1+x2+x3, data = as.data.frame(dat))
  Coef <- c(coef(f1), 0, 0)
  names(Coef)[5:6] <- c("x12", "x23")
  
  return(list(Coef=Coef, Sigma = sigma(f1)))
}

gen_pyx <- function(model_out, x, B)
{
  Coef <- model_out$Coef
  Sig <- model_out$Sigma
  
  mu <- sum(x*Coef)
  return(mu+rnorm(B, 0, Sig))
}
