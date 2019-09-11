#' Gauss-Newton update for alphas.
#'
#' Computes the update of alphas in the GAIM by in a Gauss-Newton algorithm.
#'
#' @param r Numeric vector. The residual of the current Gauss-Newton step.
#' @param x List of the index matrices. 
#' @param dgz Numeric matrix containing the derivative of current smooth
#'    functions.
#' @param alpha List of the current values of alpha coefficients.
#' @param w Numeric vector of weights.
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
#! POTENTIALLY RIDGE REGRESSION OR LASSO? SEE ROOSEN AND HASTIE (1994) AND SEARCH LITTERATURE
gn_update <- function(r, x, dgz, alpha = rep(0, sum(pvec)), 
  w = rep(1 / length(y), length(y)), delta = TRUE, monotone = 0, 
  sign.const = 0, constraint.algo = c("QP", "reparam"), norm.type = "L2",
  solver = NULL, qp_pars = list())
{
  constraint.algo <- match.arg(constraint.algo)
  if(is.matrix(x)) x <- list(x)
  p <- length(x)
  pvec <- sapply(x, ncol)
  pind <- rep(1:p, pvec)
  alphavec <- unlist(alpha)
  # prepare predictors
  dgz <- as.data.frame(dgz) 
  zerod <- apply(dgz, 2, function(x) all(x == 0))
  if (any(zerod)) dgz[,zerod] <- .Machine$double.eps
  V <- Map("*", x, dgz)
  Vmat <- data.matrix(Reduce("cbind", V))
  # Prepare constraints
  Amat <- const.matrix(pind, monotone, sign.const)
  # prepare response
  if (!delta){ 
    r <- r + Vmat %*% alphavec
    l <- rep(0, nrow(Amat)) 
  } else {
    l <- -(Amat %*% alphavec)
  }
  if (any(monotone != 0 | sign.const != 0)){   
    if (any(!monotone %in% -1:1 | !sign.const %in% -1:1)){
      stop("monotone and sign.const must be one of c(-1, 0, 1)")
    }
    alpha.up <- switch(constraint.algo,
      QP = update.QP(r, V, w, Amat, l, solver, qp_pars),
      reparam = update.reparam(r, V, w, monotone, sign.const)
    )
  } else {
    alpha.up <- coef(lm(r ~ 0 + Vmat, weights = w))
    #! To address the cases when lm cannot be fit because of important
    #!    colinearity. Naïve for now, to be thought of
    if (any(is.na(alpha.up))){
      ridge.apply <- lm.ridge(r ~ 0 + Vmat, weights = w, 
        lambda = seq(0,1,0.01))
      alpha.up <- coef(ridge.apply)[which.min(ridge.apply$GCV),]
    }
  }
  alpha.up <- split(alpha.up, pind)
#  alpha.new <- Map(normalize, alpha.up, norm.type)
  return(alpha.up)
}

#' @param x List of numeric matrices.
update.QP <- function(y, x, w, Amat, low, solver = c("osqp", "quadprog"), 
  qp_pars = list())
{
  solver <- match.arg(solver)
  p <- length(x)
  pvec <- sapply(x, ncol)
  pind <- rep(1:p, pvec)
#  pind1 <- rep(1:p, pvec - 1)
  xall <- Reduce("cbind", x)
  W <- diag(w)
  Dmat <- t(xall) %*% W %*% xall
  dvec <- 2 * t(y) %*% W %*% xall
  # OSQP
  # adaptive_rho = FALSE allows the algorithm to converge
  #   to a feasible solution.
  #   See (https://github.com/oxfordcontrol/osqp/issues/151)
  def_settings <- list(verbose = FALSE, adaptive_rho = FALSE)
  qp_pars <- c(qp_pars, 
    def_settings[!names(def_settings) %in% names(qp_pars)])
  low <- low + max(c(formals(osqp::osqpSettings)$eps_abs, qp_pars$eps_abs))
  sol <- switch(solver,
    osqp = osqp::solve_osqp(P = Dmat, q = -dvec, A = Amat, 
      l = low, u = rep(Inf, nrow(Amat)), pars = qp_pars)$x,
    quadprog = quadprog::solve.QP(Dmat, dvec, t(Amat), bvec = low)$solution
  )
  return(sol)
}

const.matrix <- function(pind, monotone, sign.const, first.pos = T){
  pind1 <- pind[c(diff(pind) == 0, FALSE)]
  Sigma <- matrix(0, length(pind) - 1, length(pind))
  diag(Sigma) <- 1
  Sigma[col(Sigma) - row(Sigma)  == 1] <- -1
  Sigma <- Sigma[diff(pind) == 0,, drop = FALSE]
  Sigma[monotone[pind1] == 1,] <- -1 * Sigma[monotone[pind1] == 1,]
  Sigma <- Sigma[monotone[pind1] != 0,, drop = FALSE]
  csign <- sign.const[pind]
  if (first.pos){
    pfirst <- which(diff(c(0, pind)) != 0)
    pfirst <- pfirst[sign.const > -1]
    csign[pfirst] <- 1
  }
  Asign <- diag(csign)
  Asign <- Asign[apply(Asign, 1, sum) != 0,, drop = F]
  return(rbind(Sigma, Asign))
}