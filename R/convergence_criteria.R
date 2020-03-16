offset_convergence <- function(r, xind, dgz)
{
  d <- ncol(xind)
  Vmat <- xind * dgz
  Vmat[abs(Vmat) < .Machine$double.eps] <- 0 # Numerical instability otherwise
  Q <- qr.Q(qr(Vmat), complete = T)  
  off <- sqrt(mean((t(Q[,1:d]) %*% r)^2) / mean((t(Q[,-(1:d)]) %*% r)^2))
  return(off)
}

L2 <- function(y, yhat, w){ 
    n <- length(y)
    if (missing(w)) w <- rep(1 / n, n)
    return(stats::weighted.mean((y-yhat)^2, w))
}
