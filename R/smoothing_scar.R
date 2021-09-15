
scar_lookup <- c(inc = "in", dec = "de", cvx = "cvx", ccv = "ccv",
  inccvx = "cvxin", deccvx = "cvxde", incccv = "ccvin",
  decccv = "ccvde"
)


scar.setup <- function(mf, control)
{
  # Vairable types
  mt <- attr(mf, "terms")
  ptot <- length(attr(mt, "term.labels"))
  gind <- attr(mt, "specials")$g
  sind <- attr(mt, "specials")$s
  ptot <- length(attr(mt, "term.labels"))
  # Create shape vector
  gsh <- rep("l", length(gind))
  fcons <- lapply(mf[gind], attr, "fcons")
  nonull <- !sapply(fcons, is.null)
  gsh[nonull] <- scar_lookup[unlist(fcons[nonull])]
  # Add covariates
  covind <- (1:ptot)[-(gind - 1)]
  csh <- rep("l", length(covind))
  fcons <- lapply(mf[covind + 1], attr, "fcons")
  nonull <- !sapply(fcons, is.null)
  csh[nonull] <- scar_lookup[unlist(fcons[nonull])]
  # Put together
  shp <- c(gsh, csh)
  if (sum(shp == "l") > 1) stop(paste0("'scar' only allows for one ",
    "unconstrained term"))
  control$shape <- shp
  return(control)
}

smooth_scar <- function(x, y, shape, Xcov = NULL, ...)
{
  scar_data <- data.matrix(cbind(x, Xcov))
  gfit <- scar::scar(x = scar_data, y = y, shape = shape, ...)
  gx <- gfit$componentfit
  dgx <- mapply(function(x, gx) stats::splinefun(x, gx)(x, deriv = 1), 
    as.data.frame(scar_data), as.data.frame(gx))
  beta0 <- gfit$constant
  return(list(intercept = beta0, gz = gx, dgz = dgx, edf = NA))
}
