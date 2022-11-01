## Compute marginal likelihood for multinomial Me
.mult_computeLogMl <- function(a, x){
  
  lbeta.xa <- sum(lgamma(a + x)) - lgamma(sum(a + x))
  lbeta.a  <- sum(lgamma(a)) - lgamma(sum(a))
  log.coef <- lgamma(sum(x) + 1) - sum(lgamma(x + 1))
  logmlHe <- (lbeta.xa-lbeta.a) + log.coef
  
  p <- 1/length(a)
  logmlH0 <- sum(x * log(p)) + log.coef
  
  output <- list(logmlHe=logmlHe,
                 logmlH0=logmlH0)
  
  return(output)
  
}

## Compute marginal likelihood for multiple binomial Me
.binom_computeLogMl <- function(a, b, x, n){
  
  lbeta.xa <- lbeta(x + a, n - x + b)
  lbeta.a  <- lbeta(a, b)
  log.coef <- lgamma(n + 1) - c(lgamma(x + 1) + lgamma(n - x + 1))
  logmlHe  <- sum(lbeta.xa - lbeta.a + log.coef)
  
  lbeta.xa.H0 <- lbeta(sum(x) + sum(a) - length(a) + 1, sum(n) - sum(x) + sum(b) - length(b) + 1)
  lbeta.a.H0  <- lbeta(sum(a) - length(a) + 1, sum(b) - length(b) + 1)
  logmlH0     <- sum(log.coef) + lbeta.xa.H0 - lbeta.a.H0 
  
  output <- list(logmlHe=logmlHe,
                 logmlH0=logmlH0)
  
  return(output)
  
}



