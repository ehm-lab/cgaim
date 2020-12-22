#' Defining indices in CGAIM
#'
#' Function used to define terms whithin a \code{cgaim} formula. \code{g} for an index on a list of variable and \code{s} for a smooth covariate.
#'
#' @param ... Variables on which the index is based. May contain vectors and matrices.
#' @param label A character used to give a name to the index. Useful for
#'    displaying the results. By default use the name of the variable passed (the first one for \code{g}).
#' @param acons A list of constraints to be applied to the index weights \code{alpha}. Can be named values for common constraints as proposed in \code{\link{build_constraints}} for built-in constraints. Can also contain numeric vectors specifying custom linear constraints, see \code{\link{alpha.setup}}. TO BE DONE
#' @param fcons The type of shape constraint to be applied on the smooth function. See details for the list of shape constraints allowed.
#' @param s_opts A named list of options to be passed to the smoothing of ridge functions. Arguments allowed heavily depend on the method used to smooth additive models. See details.
#'
#' For now, eight shape-constraints are allowed for either \code{g} or \code{s}: monotone increasing (\code{fcons = "inc"}), monotone decreasing (\code{fcons = "dec"}), convex (\code{fcons = "cvx"}), concave (\code{fcons = "ccv"}), increasing and convex(\code{fcons = "inccvx"}), decreasing and convex (\code{fcons = "deccvx"}), increasing and concave (\code{fcons = "incccv"}), decreasing and concave (\code{fcons = "decccv"}). 
#'
#' The \code{s_opts} argument can be used to pass a list of argument for basis functions used in the chosen shape-consrained smoothing method. The possible arguments when \code{smooth_method = "scam"} can be found in \code{\link[mgcv]{s}}. For \code{smooth_method = "cgam"}, the parameters allowed may vary according to the shape-contraint chosen. The full list can be found in \code{\link[cgam]{cgam}}. Only the constraints beginning with \code{s.} are allowed for now. Finally, no parameter can be passed when \code{smooth_method = "scar"} since the method does not use basis functions.
#'
#' @export
g <- function(..., label = term[1], acons = list(), 
  fcons = NULL, s_opts = list())
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
  if (!is.null(fcons)){ 
    fcons <- match.arg(fcons, c("inc", "dec", "cvx", "ccv", 
      "inccvx", "deccvx", "incccv", "decccv"))
  }
  obt_opts <- deparse(substitute(s_opts))
  opts_form <- substr(obt_opts, 6, nchar(obt_opts) - 1)
  acons <- acons[names(acons) %in% 
    methods::formalArgs(build_constraints)]
  acons$nvars <- p
  Cmat <- do.call(build_constraints, acons)
  attributes(Xmat) <- c(attributes(Xmat), 
    list(term = term, nterms = ncol(Xmat), label = label, 
    fcons = fcons, s_formula = opts_form, Cmat = Cmat))
  return(Xmat)
}

#' @describeIn g Additional smooth terms
#'
#' @param x Covariate on which the smooth is applied.
#' 
#' @export
s <- function(x, fcons = NULL, s_opts = list()){
  cl <- match.call()
  if (!is.null(fcons)){ 
    fcons <- match.arg(fcons, c("inc", "dec", "cvx", "ccv", 
      "inccvx", "deccvx", "incccv", "decccv"))
  }
  obt_opts <- deparse(substitute(s_opts))
  opts_form <- substr(obt_opts, 6, nchar(obt_opts) - 1)
  attributes(x) <- c(attributes(x),
    list(fcons = fcons, s_formula = opts_form, label = deparse(cl$x)))
  return(x)
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
  first.const = 0)
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