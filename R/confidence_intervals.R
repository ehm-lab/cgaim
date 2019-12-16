confint_gn <- function(object, parm, level = 0.95)
{
  index <- object$index
  d <- length(index)
  dp <- length(parm)
  yhat <- cbind(1, object$gz) %*% object$beta
  r <- object$y - object$fitted
  n <- length(r)
  xind <- object$x[,1:d]
  dgz <- object$dgz[,index]
  Vmat <- xind * dgz
  seV <- diag((t(Vmat) %*% Vmat)^-1)
  ssq <- sum(r^2) / (n - d)
  seAlpha <- sqrt(seV[parm]) * ssq
  tlims <- qt(c((1 - level) / 2, 1 - (1 - level) / 2), n - d)
  conflims <- sapply(seAlpha, "*", tlims)
  confs <- matrix(object$alpha, 2, dp, byrow = T) + conflims
  return(t(confs))
}

confint_bootstrap <- function(object, parm, level = 0.95)
{
# When resampled, fit the model with weighted observations.    
}