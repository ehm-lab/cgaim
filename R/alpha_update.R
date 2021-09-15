

alpha_update <- function(y, x, w, index, dgz, alpha,
  delta = TRUE, Cmat = NULL, bvec = NULL, solver, ctol, qp_pars)
{
  d <- ncol(x)
  lvec <- as.vector(table(index))
  # prepare predictors
  zerod <- apply(dgz, 2, function(x) all(x == 0))
  if (any(zerod)) dgz[,zerod] <- .Machine$double.eps
  dgz <- dgz[,index]
  Vmat <- x * dgz
  # prepare response
  if (!delta){ 
    y <- y + Vmat %*% alpha
  }
  # Normal matrices
  Dmat <- crossprod(Vmat * sqrt(w))
  dvec <- drop(crossprod(Vmat, w * y))
  # Fit
  if (!is.null(Cmat)){
    if (delta) bvec <- bvec + drop(-(Cmat %*% alpha))
    bvec <- bvec + ctol
    # bvec <- bvec + max(c(formals(osqp::osqpSettings)$eps_abs, qp_pars$eps_abs))
    res <- do.call(sprintf("update_%s", solver), 
      list(Dmat = Dmat, dvec = dvec, Cmat = Cmat, bvec = bvec, 
        qp_pars = qp_pars))
    res$active[((Cmat %*% res$alpha) - bvec) <= ctol] <- TRUE
  } else {
    # alpha.up <- stats::coef(stats::lm(y ~ 0 + Vmat, weights = w))
    alpha.up <- drop(qr.solve(Dmat) %*% dvec)
    if (any(is.na(alpha.up))){ # For non-convergence cases
      ridge.apply <- MASS::lm.ridge(y ~ 0 + Vmat, weights = w, 
        lambda = seq(0,1,0.01))
      alpha.up <- stats::coef(ridge.apply)[which.min(ridge.apply$GCV),]
    }
    res <- list(alpha = alpha.up, active = rep(FALSE, d))
  }
  return(res)
}


# update.QP <- function(y, x, w, Cmat, low, solver, qp_pars)
# {
#   Dmat <- crossprod(x * sqrt(w))
#   dvec <- drop(crossprod(x, w * y))
#   low <- low + max(c(formals(osqp::osqpSettings)$eps_abs, qp_pars$eps_abs))
#   sol <- switch(solver,
#     osqp = osqp::solve_osqp(P = Dmat, q = -dvec, A = Cmat, 
#       l = low, u = rep(Inf, nrow(Cmat)), pars = qp_pars)$x,
#     quadprog = quadprog::solve.QP(Dmat, dvec, t(Cmat), bvec = low)$solution
#   )
#   return(sol)
# }


normalize <- function(alpha, type){
  anorm <- norm(matrix(alpha, ncol = 1), type)
  if (!(is.na(anorm) || all(alpha == 0))) alpha <- alpha / anorm
  attr(alpha, "norm") <- anorm
  return(alpha)
}


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

