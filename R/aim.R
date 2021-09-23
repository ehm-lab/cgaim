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
  result <- do.call(gaim_gn, algo_pars)
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

gaim_gn <- function(x, y, w, index, 
  smooth_control = list(), alpha_control = list(),
  trace = FALSE, max.iter = 50, tol = 1e-3, min.step.len = 0.1, 
  halving = T, convergence_criterion = "rss")
{
  # Useful objects
  n <- length(y)
  p <- max(index)
  d <- length(index)
  ind_pos <- split(1:d, index)
  # Indice computation
  alpha <- alpha_control$alpha
  zs <- sapply(ind_pos, function(i, x, a) x[,i] %*% a[i], 
    x = x, a = alpha)
  colnames(zs) <- unique(names(index))
  # Initial smoothing
  smooth_fun <- sprintf("smooth_%s", smooth_control$method)
  smooth_control$method <- NULL
  smo_par <- c(smooth_control, 
    list(y = y, x = zs, weights = w))
  gz <- do.call(smooth_fun, smo_par)   
  yhat <- gz$intercept + rowSums(gz$gz)
  r <- y - yhat  
  # Convergence criterion
  eps <- tol + 1
  c1 <- 1      
  l2 <- L2(y, yhat, w)  
  # If requested: tracing the algorithm evolution
  if (trace){
    trace.list <- list(
      criterion = rep(NA, max.iter),
      alpha = matrix(NA, nrow = max.iter, ncol = d),
      gfit = rep(list(NA), max.iter),
      step.len = rep(NA, max.iter) 
    )
    trace.list$criterion[1] <- eps
    trace.list$step.len[1] <- 1  
    trace.list$alpha[1,] <- alpha 
    trace.list$gfit[[1]] <- gz$gz
  }
  # Gauss-Newton search
  alpha_pars <- c(
    alpha_control[names(alpha_control) %in% methods::formalArgs(alpha_update)],
    list(x = x, w = w, index = index, delta = TRUE))
  stopflag <- 0
  while(stopflag == 0){
    # Update alphas
    alpha_pars$y <- r
    alpha_pars$dgz <- gz$dgz
    alpha_pars$alpha <- alpha
    alphaup <- do.call(alpha_update, alpha_pars)
    delta <- alphaup$alpha
    # Halving in case of bad steps
    cuts <- 1
    repeat {   
      alpha.new <- alpha + delta * cuts
      alpha.new <- unlist(tapply(alpha.new, index, normalize, 
        type = alpha_control$norm.type))
      zs <- sapply(ind_pos, function(i, x, a) x[,i] %*% a[i], 
        x = x, a = alpha.new)
      colnames(zs) <- unique(names(index))   
      smo_par$x <- zs
      gz <- do.call(smooth_fun, smo_par)
      yhat <- gz$intercept + rowSums(gz$gz)
      l2.new <- L2(y, yhat, w)
      if (l2 - l2.new > 0 || !halving){
         break
      }
      cuts <- cuts / 2
      if (cuts < min.step.len){ 
        stopflag <- 2
        break
      }
    }
    if (convergence_criterion == "offset"){
      eps <- offset_convergence(r, x, gz$dgz[,index])
    } else if (convergence_criterion == "alpha") {
      eps <- max(abs(alpha.new - alpha) / abs(alpha))
    } else {
      eps <- (l2 - l2.new) / l2
    }
    r <- y - yhat
    alpha <- alpha.new
    l2 <- l2.new
    c1 <- c1 + 1       
    if (trace){ # Tracing the algorithm evolution
      trace.list$criterion[c1] <- eps
      trace.list$step.len[c1] <- cuts 
      trace.list$alpha[c1,] <- alpha 
      trace.list$gfit[[c1]] <- gz$gz
    }
    if (all(eps < tol)){
      stopflag <- 1
    } else { 
      if (c1 >= max.iter){
        stopflag <- 3
        warning(paste0("Fitting did not converge after ", max.iter, 
          " iterations. Consider revising the constraints"))
      } 
    }
  }
  names(alpha) <- colnames(x)
  final.gz <- scale(gz$gz)
  final.gz[,apply(gz$gz, 2, stats::sd) == 0] <- 
    gz$gz[,apply(gz$gz, 2, stats::sd) == 0]
  betas <- attr(final.gz, "scaled:scale")
  r <- y - yhat
  edf <- gz$edf + d - sum(alphaup$active)
  gcv <- l2 / (1 - (edf / n))^2
  sig2 <- sum(r^2) / (n - edf)
  output <- list(alpha = alpha, gfit = final.gz, indexfit = zs, 
    beta = c(gz$intercept + sum(attr(final.gz, "scaled:center")), betas),
    index = index, fitted = yhat, residuals = r, rss = l2, flag = stopflag, 
    niter = c1, edf = edf, gcv = gcv)
  if (trace){
    trace.list$criterion <- trace.list$criterion[1:c1]
    trace.list$alpha <- trace.list$alpha[1:c1,]
    trace.list$gfit <- trace.list$gfit[1:c1]
    trace.list$step.len <- trace.list$step.len[1:c1]
    output$trace <- trace.list
  }
  return(output)   
}