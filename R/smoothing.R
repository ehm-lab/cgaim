#' Estimate ridge functions g
#'
#' Performs scatterplot smoothing to obtain ridge functions g.
#'
#' @param x A matrix or data frame giving indices (one per column).
#' @param y A numeric vector containing the output of the model.
#' @param w A numeric vector containing weights.
#' @param method A character value giving the method for smoothing.
#'    'mgcv' calls either 'gam' or 'scam' depending on whether there are
#'    shape-constrained terms. 'scar' call 'scar'.
#' @param shape A character vector indicating the type of smoothing for each
#'    index. Recycled if not the same length as \code{ncol(x)}. Can be one of
#'    the smoothers available in the \code{mgcv} package (see 
#'    code{\link[mgcv]{smooth.terms}}). Can also be one of the 
#'    shaped-constrained smoothers in \code{scam} (see 
#'    code{\link[scam]{shape.constrained.smooth.termss}}).
#'    Can also be "l" for linear.
#' @param ... Additional arguments to be passed to the method.
#! Parameter organization in the main function (especially for shapes)
smoothing <- function(x, y, w, method = c("mgcv", "scar"), 
  shape = "tp", ...)
{
  dots <- list(...)
  method <- match.arg(method)
  p <- ncol(x)
  x <- as.data.frame(x)
  names(x) <- sprintf("V%i", 1:p)
  n <- length(y)
  shape <- rep_len(shape, p)
  smooth_pars <- c(dots, list(x = x, y = y, w = w, shape = shape))
  out <- do.call(sprintf("smooth_%s", method), smooth_pars)
  return(out)  
}

smooth_mgcv <- function(x, y, w, shape, ...)
{
  dots <- list(...)
  lin <- shape == "l"
  nlin <- sum(lin)
  form.rhs <- colnames(x)
  form.rhs[!lin] <- sprintf("s(%s, bs = '%s')", 
    form.rhs[!lin], shape[!lin])
  form <- sprintf("y ~ %s", paste(form.rhs, collapse = " + "))  
  # Determine model to fit
  mod <- ifelse(any(shape %in% 
    c("mpi", "mpd", "cx", "cv", "micx", "micv", "mdcx", "mdcv")),
    scam::scam, mgcv::gam
  )  
  mgcv_pars <- c(list(formula = as.formula(form), 
    data = data.frame(y = y, x), weights = w), 
    dots[names(dots) %in% names(formals(mod))]
  )  
  gfit <- do.call(mod, mgcv_pars)
  # Extract estimated terms
  gx <- predict(gfit, type = "terms")
  # Estimate first derivative (different functions for scam or gam)
  if (any(shape %in% 
    c("mpi", "mpd", "cx", "cv", "micx", "micv", "mdcx", "mdcv")))
  {
    dmod <- lapply(1:sum(!lin), scam::derivative.scam, object = gfit)    
    dsm <- sapply(dmod, "[[", "d")    
  } else {
    dmod <- gratia::fderiv(gfit, newdata = x)
    dsm <- sapply(dmod$derivatives, "[[", "deriv")    
  }  
  ny <- length(y)
  dgx <- matrix(NA, ny, ncol(x))
  if (sum(!lin) > 0) dgx[,!lin] <- dsm
  suppressWarnings(dgx[,lin] <- matrix(coef(gfit)[1:nlin + 1], 
      nrow = ny, ncol = nlin, byrow = TRUE))
  beta0 <- coef(gfit)[1]
  return(list(intercept = beta0, gz = gx, dgz = dgx, fit = gfit))
}

smooth_scar <- function(x, y, w, shape, ...)
{
  dots <- list(...)
  scar.pars <- c(list(x = data.matrix(x), y = y, shape = shape, weights = w), 
    dots[names(dots) %in% names(formals(scar::scar))])
  gfit <- do.call(scar::scar, scar.pars)
  gx <- gfit$componentfit
  dgx <- mapply(function(x, gx) splinefun(x, gx)(x, deriv = 1), 
    as.data.frame(x), as.data.frame(gx))
  beta0 <- gfit$constant
  return(list(intercept = beta0, gz = gx, dgz = dgx, fit = gfit))
}

