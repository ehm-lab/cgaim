#!############################################################################
#!
#!           A global AIM
#!
#!############################################################################

#' Constrained groupwise additive index models
#'
#' Fits a constrained groupwise additive index model (CGAIM) through alternative 
#' sequential quadratic programming steps.
#' 
#' The CGAIM is expressed 
#'  \deqn{y_{i} = \beta_{0} + \sum_{j} \beta_{j} g_{j}(\alpha_{j}^{T} x_{ij})
#'    + \sum_{k} \gamma_{k} f_{k}(x_{ik}) + e_{i}}
#'  The formula interface considers \code{\link{g}} to identify index terms,
#'   \code{\link{s}} for smooth functions and can also include 
#'    linear terms. All smooth terms can be shape constrained.
#'
#' The CGAIM allows for linear constraints on the alpha coefficients. 
#'  Common constraints can be given for each index through the \code{acons} 
#'  parameter in \code{\link{g}} terms in the formula. Custom constraints 
#'  can either be given for each idnex individually through the 
#'  \code{\link{g}} interface, or for all parameters through
#'  \code{alpha.control$Cmat}. See \code{\link{alpha.control}}.
#'
#' The CGAIM is fitted through an iterative algorithm that alternates between
#'  estimating the ridge functions \eqn{g_{j}} (and other non-index terms) and 
#'  updating the coefficients \eqn{\alpha_{j}}. The smoothing of ridge functions
#'  currently supports three methods: \code{\link[scam]{scam}} (the default), 
#'  \code{\link[cgam]{cgam}} and \code{\link[scar]{scar}}. 
#'  The list \code{smooth.control} allows controlling the functions.
#'
#' Updating the coefficient is made through quadratic programming. Currently,
#'  it is carried out by either the function 
#'  \code{\link[osqp:solve_osqp]{solve_osqp}} or 
#'  \code{\link[quadprog:solve.QP]{solve.QP}}. See \code{\link{alpha.control}}
#'  for details.
#'
#' @param formula A CGAIM formula with index terms \code{\link{g}}, 
#'  smooth terms \code{\link{s}} and linear terms. 
#' @param data A data.frame containing the variables of the model.    
#' @param weights A vector of optional weights for observations. 
#' @param na.action A function indicating what treatment applying to NAs. If
#'  missing the default is set by the \code{na.action} setting of \code{options}. 
#'  See \code{\link[stats]{na.fail}}.
#' @param smooth_method The shape constrained smoothing method to consider. 
#'  The default is \code{\link[scam]{scam}}, but other options are 
#'  \code{\link[cgam]{cgam}} and \code{\link[scar]{scar}}. See details.
#' @param smooth_control A list of control parameters for the shape-constrained 
#'  smoothing function given in \code{smooth_method}. 
#'  See help of the functions for the list of available parameters.
#' @param alpha_control A list of parameters controlling the
#'    alpha updating steps of the algorithm. See \code{\link{alpha.control}}.
#' @param algo_control A list containing controlling parameters for the
#'    outer algorithm. See \code{\link{algo.control}}.
#'
#' @return A \code{cgaim} object, i.e. a list with components:
#'  \item{alpha}{A named list of index coefficients.}
#'  \item{gfit}{A matrix containing the ridge and smooth functions 
#'    evaluated at the observations. 
#'    Note that column ordering puts indices first and covariates after.}
#'  \item{dg}{A matrix containing derivatives of ridge and smooth functions.}
#'  \item{indexfit}{A matrix containing the indices evaluated at the 
#'    observations.}
#'  \item{beta}{A vector containing the intercept and the scale coefficient
#'    of each ridge and smooth function. Includes the \eqn{\gamma_{k}} of
#'    the CGAIM model above. Note that ordering puts indices first and 
#'    covariates after.}
#'  \item{index}{A vector identifying to which index the columns of the
#'    element \code{x} belong.}
#'  \item{fitted}{A vector of fitted responses.}
#'  \item{residuals}{A vector of residuals.}
#'  \item{rss}{The residual sum of squares of the fit.}
#'  \item{flag}{A flag indicating how the algorithm stopped. 1 for proper 
#'    convergence, 2 when the algorithm stopped for failing to decrease rss
#'    and 3 when the maximum number of iterations has
#'    been reached.}
#'  \item{niter}{Number of iterations performed.}
#'  \item{edf}{Effective degrees of freedom of the estimator.}
#'  \item{gcv}{Generalized cross validation score.}
#'  \item{covariates}{A data.frame containing the values of the variables not 
#'    entering any index.}
#'  \item{x}{A matrix containing the variables entering the indices. 
#'    The variables are mapped to each index through the element \code{index}.}
#'  \item{y}{The response vector.}
#'  \item{weights}{The weights used for estimation.}
#'  \item{smooth_control}{List of parameters controlling the smoothing step.}
#'  \item{alpha_control}{List of parameters controlling alpha updating step.}
#'  \item{algo_control}{List of parameters controlling the algorithm.}
#'  \item{terms}{A \code{\link[stats]{terms.object}} representing the model.
#'    Useful for prediction.}
#'  \item{trace}{If \code{algo_control$trace = TRUE}, a tracking list 
#'    containing, for each iteration, values of alpha coefficients, 
#'    ridge and smooth function, convergence criteria, step length,
#'    as well as the number of iterations.}
#'
#' @seealso \code{\link{confint.cgaim}} for confidence interval,
#'    \code{\link{predict.cgaim}} to predict new data,
#'    \code{\link{plot.cgaim}} to plot ridge function.
#'
#' @export
cgaim <- function(formula, data, weights, na.action, 
  smooth_method = "scam", smooth_control = list(), 
  alpha_control = list(), algo_control = list()) 
{
  # Terms
  mt <- stats::terms(formula, specials = c("g", "s"), data = data)
  gind <- attr(mt, "specials")$g
  p <- length(gind)
  ptot <- length(attr(mt, "term.labels"))
  # Fill missing arguments
  # n <- NROW(data[[1]])
  # if (is.null(weights)){
  #   weights <- rep(1, n)
  # } else {
  #   weights <- rep_len(weights, n)
  # }
  # if (is.null(na.action)) na.action <- getOption("na.action")
  # Extract data
  cl <- match.call(expand.dots = FALSE)
  m <- match(c("data", "weights", "na.action"), 
    names(cl), 0L)
  cl <- cl[c(1L, m)]
  cl$drop.unused.levels <- TRUE
  cl$formula <- mt
  cl[[1L]] <- quote(stats::model.frame)
  mf <- eval(cl, parent.frame())
  cl$formula <- stats::reformulate(all.vars(mt))
  refordata <- eval(cl, parent.frame())
  # mf <- stats::model.frame(mt, data = data, 
  #   na.action = na.action, weights = weights)
  # Initialize index estimation components
  mod_components <- index.setup(mf, alpha_control)
  # Prepare parameters for smoothing
  smooth_components <- smooth.setup(mf, refordata, 
    smooth_method, smooth_control)
  # Prepare parameters for fitting algorithm
  algo_control <- do.call(algo.control, algo_control)
  # Fitting the model
  algo_pars <- c(mod_components, algo_control, 
    list(smooth_control = smooth_components))
  result <- do.call(cgaim.fit, algo_pars)
  # Organize output
  result$alpha <- split(result$alpha, mod_components$index)
  names(result$alpha) <- unique(names(mod_components$index))
  attributes(result$gfit) <- attributes(result$gfit)[c("dim", "dimnames")]
  result$covariates <- smooth_components$Xcov
  result$x <- mod_components$x
  result$y <- mod_components$y
  result$weights <- mod_components$w
  result$smooth_control <- smooth_components
  result$alpha_control <- mod_components$alpha_control
  result$algo_control <- algo_control
  result$terms <- mt
  class(result) <- "cgaim"
  return(result)
}

