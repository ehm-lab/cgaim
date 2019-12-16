#' @param type Character string indicating the type of initialization.
#'    'regression' regresses x on y, 'random' draws alphas from a uniform
#'    distribution and 'constant' just set them at 1 / p

alpha_init <- function(y, x, w, index, 
  init.type = c("regression", "random"), Cmat = NULL, 
  norm.type = "L2", ...)
{
  init.type <- match.arg(init.type)
  dots <- list(...)
  alpha <- switch(init.type,
    random = {
      pars <- dots[names(dots) %in% names(formals(limSolve::xsample))]
      pars$G <- Cmat
      pars$H <- rep(0, nrow(Cmat))
      res <- suppressWarnings(do.call(limSolve::xsample, pars))$X
      res[nrow(res),]
    },
    regression = {
      pars <- dots[names(dots) %in% names(formals(alpha_update))]
      pars <- within(pars, {
        y <- y; x <- x; w <- w; index <- index; Cmat <- Cmat
        dgz <- matrix(1, length(y), max(index))
        alpha <- rep(0, ncol(x))
        delta <- FALSE
      })
      do.call(alpha_update, pars)
    }
  )
  # Normalization
  alpha <- unlist(tapply(alpha, index, normalize, type = norm.type))
  return(alpha)
}


#' Gauss-Newton update for alphas.
#'
#' Computes the update of alphas in the GAIM by in a Gauss-Newton algorithm.
#'
#' @param y Numeric vector. The response of the current step.
#' @param x Matrix of all variables.
#' @param w Numeric vector of weights.
#' @param index The indices 
#' @param dgz Numeric matrix containing the derivative of current smooth
#'    functions.
#' @param alpha List of the current values of alpha coefficients.
#' @param delta Logical. If TRUE, the change is computed instead of the
#'    new alphas directly.
#' @param monotone Integer vector indicating if alphas are constrained to be
#'    monotone. 0 indicates no monotonicity, -1 decreasing and 1 increasing.
#'    Recycled if of a different length than dg.
#' @param sign.const Integer vector indicating a constraint on the sign of
#'    alphas. 0 indicates no constraint, -1 that all alphas must be negative 
#'    and 1 that all alphas must be positive. Recycled if of a different length 
#'    than dg.
#' @param constraint.algo Character indicating which method is used to contraint 
#'    monotonicity of parameters. "QP" indicates classical quadratic programming
#'    and "pya" indicates the method used in the SCAM (Pya and Wood, 2015).
alpha_update <- function(y, x, w, index, dgz, alpha, 
  delta = TRUE, Cmat = NULL, solver = NULL, qp_pars = list())
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
  if (!is.null(Cmat)){
    l <- if (delta) drop(-(Cmat %*% alpha)) else rep(0, nrow(Cmat))
    alpha.up <- update.QP(y, Vmat, w, Cmat, l, solver, qp_pars)    
  } else {
    alpha.up <- coef(lm(y ~ 0 + Vmat, weights = w))
    if (any(is.na(alpha.up))){ # For non-convergence cases
      ridge.apply <- lm.ridge(y ~ 0 + Vmat, weights = w, 
        lambda = seq(0,1,0.01))
      alpha.up <- coef(ridge.apply)[which.min(ridge.apply$GCV),]
    }
  }
  return(alpha.up)
}

#' @param x List of numeric matrices.
update.QP <- function(y, x, w, Cmat, low, solver = c("osqp", "quadprog"), 
  qp_pars = list())
{
  solver <- match.arg(solver) 
  W <- diag(w)
  Dmat <- t(x) %*% W %*% x
  dvec <- 2 * t(y) %*% W %*% x
  # OSQP
  # adaptive_rho = FALSE allows the algorithm to converge
  #   to a feasible solution.
  #   See (https://github.com/oxfordcontrol/osqp/issues/151)
  def_settings <- list(verbose = FALSE, adaptive_rho = FALSE)
  qp_pars <- c(qp_pars, 
    def_settings[!names(def_settings) %in% names(qp_pars)])
  low <- low + max(c(formals(osqp::osqpSettings)$eps_abs, qp_pars$eps_abs))
  sol <- switch(solver,
    osqp = osqp::solve_osqp(P = Dmat, q = -dvec, A = Cmat, 
      l = low, u = rep(Inf, nrow(Cmat)), pars = qp_pars)$x,
    quadprog = quadprog::solve.QP(Dmat, dvec, t(Cmat), bvec = low)$solution
  )
  return(sol)
}

const.matrix <- function(index, monotone, sign.const, first.pos = T){
  index1 <- index[c(diff(index) == 0, FALSE)]
  Sigma <- matrix(0, length(index) - 1, length(index))
  diag(Sigma) <- 1
  Sigma[col(Sigma) - row(Sigma)  == 1] <- -1
  Sigma <- Sigma[diff(index) == 0,, drop = FALSE]
  Sigma[monotone[index1] == 1,] <- -1 * Sigma[monotone[index1] == 1,]
  Sigma <- Sigma[monotone[index1] != 0,, drop = FALSE]
  csign <- sign.const[index]
  if (first.pos){
    pfirst <- which(diff(c(0, index)) != 0)
    pfirst <- pfirst[sign.const > -1]
    csign[pfirst] <- 1
  }
  Asign <- diag(csign)
  Asign <- Asign[apply(Asign, 1, sum) != 0,, drop = F]
  return(rbind(Sigma, Asign))
}

#' Normalization of a vector of alphas
#'
#' Provides different ways to normalize a vector.
#'
#' @param alpha Numeric vector to normalize.
#' @param type Character indicating the type of normalization. List...
normalize <- function(alpha, type = "L2"){
  anorm <- switch(type,
    L2 = norm(alpha, "2"),
    L1 = norm(as.matrix(alpha), "1"),
    sum = sum(alpha),
    1
  )
  out <- if (anorm != 0) alpha / anorm else alpha
  attr(out, "norm") <- anorm
  return(out)
}