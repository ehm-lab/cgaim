#' Defining indices in CGAIM
#'
#' Function used to define indices whithin a \code{cgaim} formula. 
#'
#' @param ... A list of variables on which the index is based. May contain vectors and matrices.
#' @param label A character used to give a name to the index. Useful for
#'    displaying the results. 
#' @param constraints A list of common constraints for the index. See 
#'    \code{\link{build_constraints}} for built-in constraints. For
#'    other constraints, see \code{\link{alpha.setup}}.
#' @param bs The type of basis functions for the ridge function. Used to
#'    give shape constraints. The bases allowed are expressed as in the
#'    \code{scam} (\code{\link[scam]{shape.constrained.smooth.terms}}) and
#'    \code{mgcv} (\code{\link[mgcv]{smooth.terms}}) packages.
#' @param sOptions A named list of options to be passed to \code{\link[mgcv]{s}}
#'    for the ridge function smoothing.
#'
#' @export
g <- function(..., label = term[1], constraints = list(), bs = "tp", 
  sOptions = list())
{
  vars <- as.list(substitute(list(...)))[-1]
  term <- sapply(vars, deparse)
  Xvars <- list(...)
  for (j in seq_len(length(Xvars))){
    Xvars[[j]] <- as.matrix(Xvars[[j]])
    if (ncol(Xvars[[j]]) > 1 && is.null(colnames(Xvars[[j]]))) 
      colnames(Xvars[[j]]) <- as.character(seq_len(ncol(Xvars[[j]])))
  }
  pvec <- sapply(Xvars, ncol)
  n <- unique(sapply(Xvars, nrow))
  if (length(n) > 1) stop("All variables in 'g' must have the same length")  
  Xmat <- do.call(cbind, Xvars)
  colnames(Xmat) <- unlist(Map(paste0, term, lapply(Xvars, colnames)))
  p <- ncol(Xmat)
  if (p < 2) warning(sprintf("Less than 2 variables for index %s", label))
  s_formula <- sprintf("bs = \"%s\"", bs)
  obt_opts <- deparse(substitute(sOptions))
  opts_form <- substr(obt_opts, 6, nchar(obt_opts) - 1)
  if (nchar(opts_form) > 0) 
    s_formula <- paste(c(s_formula, opts_form), collapse = ", ")
  constraints <- constraints[names(constraints) %in% 
    methods::formalArgs(build_constraints)]
  constraints$nvars <- p
  Cmat <- do.call(build_constraints, constraints)
  attributes(Xmat) <- c(attributes(Xmat), 
    list(term = term, nterms = ncol(Xmat), label = label, 
    opt_formula = s_formula, Cmat = Cmat))
  return(Xmat)
}

#' Common constraints
#'
#' Build a constraint matrix from common simple constraints.
#' 
#' @param nvars A vector giving, for each index, its number of variables.
#' @param monotone A vector the same length as \code{nvars} giving
#'    monotonicity constraints. \code{0} means no constraint, \code{-1}
#'    decreasing coefficients and \code{1} increasing coefficients.
#' @param sign.const A vector the same length as \code{nvars} giving
#'    signs constraints. \code{0} means no constraint, \code{-1}
#'    means that, for the corresponding index, all \eqn{\alpha <= 0} and
#'    \code{-1} all \eqn{\alpha >= 0}.
#' @param first.const A vector indicating a sign constraint for first
#'    coefficient of each index. Recommended for indentifiability purposes but
#'    overriden by \code{sign.const}.
#'
#' @export
build_constraints <- function(nvars, monotone = 0, sign.const = 0, 
  first.const = 1)
{
  p <- length(nvars)
  monotone <- rep_len(monotone, p)
  sign.const <- rep_len(sign.const, p)
  first.const <- rep_len(first.const, p)
  index <- rep(1:p, nvars)
  index1 <- rep(1:p, nvars - 1)
  Sigma <- matrix(0, length(index) - 1, length(index))
  diag(Sigma) <- 1
  Sigma[col(Sigma) - row(Sigma)  == 1] <- -1
  Sigma <- Sigma[diff(index) == 0,, drop = FALSE]
  Sigma[monotone[index1] == 1,] <- -1 * Sigma[monotone[index1] == 1,]
  Sigma <- Sigma[monotone[index1] != 0,, drop = FALSE]
  csign <- sign.const[index]
  pfirst <- which(diff(c(0, index)) != 0)
  first.const[sign.const != 0] <- sign.const[sign.const != 0]
  csign[pfirst] <- first.const
  Asign <- diag(csign)
  Asign <- Asign[apply(Asign, 1, sum) != 0,, drop = F]
  return(rbind(Sigma, Asign))
}


#' Iterative algorithm control
#'
#' Internal function setting up the control of the estimation algorithm.
#'  Use this as a comprehensive list of allowed parameters for the 
#'  \code{algo.control} argument in \code{\link{cgaim}}.
#'
#' @param max.iter The maximum number of iteration allowed. Default to 50.
#' @param tol The tolerance for convergence. The algorithm stops when 
#'    the convergence criteria fall below \code{tol}.
#' @param min.step.len The minimum descent step length to be considered. Each
#'    step has a default step length of 1, but it is halved if it results
#'    in a step such that the convergence criterion increases. It is halved
#'    until resulting in a descent step of falling below 
#'    \code{min.step.len}.
#' @param halving Logical indicating if halving the step length in case of bad 
#'    step should be performed.
#' @param convergence_criterion The criterion used to track convergence of the
#'    algorithm. When \code{convergence_criterion = "change"}, the algorithm
#'    stops when the residual sum of squares and alpha coefficients change
#'    by less than \code{tol}. Note that, in this case \code{tol} can be
#'    a vector with two values. When \code{convergence_criterion = "offset"},
#'    the offset criterion is used. It measures the orthogonality between 
#'    the descent direction and residual sum of squares. This is the 
#'    criterion used for instance by the \code{link[stats]{nls}} function.
#'
#' @references
#'    Bates, D.M., Watts, D.G., 1981. A Relative Off set 
#'      Orthogonality Convergence Criterion for Nonlinear least Squares. 
#'      Technometrics 23, 179–183.
#'
#' @export
algo.setup <- function(max.iter = 50, tol = 1e-3, min.step.len = 0.1, 
  halving = T, convergence_criterion = c("change", "offset"))
{
  invisible(NULL)  
}