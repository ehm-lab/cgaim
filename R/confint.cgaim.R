#' Confidence intervals
#'
#' Computes confidence intervals for the CGAIM parameters.
#'    
#' @param object A \code{cgaim} object.
#' @param parm The type of parameters for which to get confidence intervals.
#' @param level The level of confidence intervals. Default to 95\%.
#' @param type A character vector giving the type(s) of confidence intervals
#'    to compute. \code{"normal"} compute classical confidence intervals
#'    based on normal approximation. The others are all based on bootstrap.
#'    \code{"boot.pct"} computes intervals as percentiles of the bootstrap
#'    distribution. \code{"boot.t"} computes studentized bootstrap intervals.
#'    \code{"boot.bca"} (still experimental, use at your own risk) 
#'    computes bias-corrected and accelerated bootstrap
#'    intervals. Can include several types.
#'    CURRENTLY, ONLY "boot.pct" WORKS
#' @param ... Additional parameters to be passed to the parallel function. See
#'    \code{\link[parallel]{mclapply}} or \code{\link[parallel]{parLapply}}.
#'
#' @return A named list with elements \code{alpha}, \code{beta} and \code{g}.
#'    Each is a list containing as elements the different types of 
#'    intervals given in \code{type}. In addition, contains the generated
#'    bootstrap values as well as the bias and acceleration computed if
#'    \code{"boot.bca"} is in \code{type}.
#'
#' @references
#'   Pya, N., Wood, S.N., 2015. Shape constrained additive models. 
#'    Stat. Comput. 25, 543–559. 
#'   Wood, S.N., 2017. Generalized Additive Models: An Introduction with R, 
#'     2nd ed, Texts in Statistical Science. Chapman and Hall/CRC.
#'   DiCiccio, T.J., Efron, B., 1996. Bootstrap Confidence Intervals. 
#'     Statistical Science 11, 189–212.
#'   Carpenter, J., Bithell, J., 2000. Bootstrap confidence intervals: 
#'     when, which, what? A practical guide for medical statisticians. 
#'     Statistics in Medicine 19, 1141–1164.
#'
#' @export
confint.cgaim <- function(object, parm, level = 0.95, 
  type = c("normal", "bootstrap"), B = 100, ...)
{
  #----- header
  
  # Check type
  type <- match.arg(type)
  
  # Get useful objects
  p <- length(object$alpha)
  ptot <- ncol(object$gfit)
  d <- length(object$index)
  n <- length(object$fitted)
  
  # Setup parm. If missing get all parameters
  if (missing(parm)){
    parm <- 1:3
  } else {      
    if (is.character(parm)){
      parm <- na.omit(match(parm, c("alpha", "beta", "g")))
    } else {
      parm <- parm[parm %in% 1:3]
    }
    if (!any(parm %in% 1:3)){
      stop(paste0("'parm' must be in 1:3"))
    }
  }
  
  # Setup CI level
  alims <- c((1 - level) / 2, 1 - (1 - level) / 2)
  level.labels <- sprintf("%s%%", alims * 100)
  
  #----- Compute confidence intervals
  if (type == "normal"){
    
    res <- list()
    
    # Alpha
    if (1 %in% parm){
      # Simulate alpha from truncated multivariate normal
      simures <- simul_tmvnorm(object, B = B)
      
      # Compute CI for Alpha
      res$alpha <- t(apply(simures, 2, stats::quantile, alims, na.rm = TRUE))
      
      # Names
      rownames(res$alpha) <- names(unlist(object$alpha))
      colnames(res$alpha) <- level.labels
    }
    
    # Beta
    if (2 %in% parm){
      # Extract standard error for beta
      vb <- vcov_beta(object)
      seb <- sqrt(diag(vb))
      
      # Compute CI for beta
      betahat <- object$beta
      res$beta <- betahat + matrix(seb, ncol = 1) %*% qnorm(alims)
      
      # names
      rownames(res$beta) <- names(object$beta)
      colnames(res$beta) <- level.labels
    }
    
    # Ridge functions
    if (3 %in% parm){
      # Extract bounds
      gbnd <- qnorm(alims[2]) * t(t(object$gse) / object$beta[-1])
      
      # Boundaries
      res$g <- array(NA, dim = c(n, ptot, 2), 
        dimnames = list(NULL, colnames(object$gfit), level.labels))
      res$g[,,1] <- object$gfit - gbnd
      res$g[,,2] <- object$gfit + gbnd
    }
    
  } else {
    
    # Simulate
    simures <- boot.cgaim(object, B = B, ...)
    
    # Compute CI
    res <- confint(simures, parm = parm, level = level)
  }
  
  #----- Return
  return(res)
}
