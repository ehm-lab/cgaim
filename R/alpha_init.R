
alpha_init <- function(y, x, w, index, 
  init.type = c("regression", "random"), Cmat = NULL, bvec = NULL,
  norm.type = "1", ...)
{
  init.type <- match.arg(init.type)
  dots <- list(...)
  alpha <- switch(init.type,
    random = {
      pars <- dots[names(dots) %in% names(formals(limSolve::xsample))]
      pars$G <- Cmat
      pars$H <- bvec
      res <- suppressWarnings(do.call(limSolve::xsample, pars))$X
      res[nrow(res),]
    },
    regression = {
      pars <- dots[names(dots) %in% names(formals(alpha_update))]
      pars <- within(pars, {
        y <- y; x <- x; w <- w; index <- index; Cmat <- Cmat; bvec <- bvec
        dgz <- matrix(1, length(y), max(index))
        alpha <- rep(0, ncol(x))
        delta <- FALSE
      })
      do.call(alpha_update, pars)$alpha
    }
  )
  # Normalization
  alpha <- unlist(tapply(alpha, index, normalize, type = norm.type))
  return(alpha)
}