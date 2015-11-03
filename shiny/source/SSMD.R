ssmd <- function(pop1, pop2, alpha=0.05)
{
  fctr <- function(pop, exp.denominator=1, exp.sd=2)
  {
    return ((length(pop) - 1) / (length(pop))^exp.denominator * (sd(pop))^exp.sd)
  }
  beta <- mean(pop1) - mean(pop2)
  beta <- beta / sqrt(fctr(pop1) + fctr(pop2))
  
  errorA <- ( fctr(pop1, exp.denominator=2) + fctr(pop2, exp.denominator=2) ) / ( fctr(pop1) + fctr(pop2) )
  errorB <- fctr(pop1, exp.denominator=3, exp.sd=4) + fctr(pop2, exp.denominator=3, exp.sd=4)
  errorC <- ( mean(pop1) - mean(pop2) )^2 / (2 * (fctr(pop1) + fctr(pop2))^3)
  
  error <- pnorm(alpha) * sqrt(errorA + (errorB * errorC))
  
  return (c(beta, error))
}

ssmd.robust <- function(pop1, pop2, alpha=0.05)
{
  fctr <- function(pop, exp.denominator=1, exp.mad=2)
  {
    return ((length(pop) - 1) / (length(pop))^exp.denominator * (mad(pop))^exp.mad)
  }
  beta <- median(pop1) - median(pop2)
  beta <- beta / sqrt(fctr(pop1) + fctr(pop2))
  
  errorA <- ( fctr(pop1, exp.denominator=2) + fctr(pop2, exp.denominator=2) ) / ( fctr(pop1) + fctr(pop2) )
  errorB <- fctr(pop1, exp.denominator=3, exp.mad=4) + fctr(pop2, exp.denominator=3, exp.mad=4)
  errorC <- ( median(pop1) - median(pop2) )^2 / (2 * (fctr(pop1) + fctr(pop2))^3)
  
  error <- pnorm(alpha) * sqrt(errorA + (errorB * errorC))
  
  return (c(beta, error))
}

ssmdMM <- function(pop1, pop2)
{
  beta <- mean(pop1) - mean(pop2)
  beta <- beta / sqrt((sd(pop1))^2 + (sd(pop2))^2)
  return (beta)  
}