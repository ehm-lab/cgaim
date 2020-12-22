#!############################################################################
#!
#!           A global AIM
#!
#!############################################################################

#' Constrained groupwise additive index models
#'
#' Fits a constrained groupwise additive index model (CGAIM) through
#'  alternative sequential quadratic programming steps.
#'
#' @param formula A CGAIM formula with index terms \code{\link{g}}, 
#'    smooth terms \code{\link{s}} and linear terms. 
#' @param data A data.frame containing the variables of the model.    
#' @param weights A vector of optional weights for observations. If \code{NULL}, unit weights are considered.
#' @param na.action A function indicating what treatment applying to NAs. If
#'    missing the default is set by the \code{na.action} setting of \code{options}. See
#'    \code{\link[stats]{na.fail}}.
#' @param smooth_method The shape constrained smoothing method to consider. The default is \code{\link[scam]{scam}}, but other options are \code{\link[cgam]{cgam}} and \code{\link[scar]{scar}}. See details.
#' @param smooth.control A list containing parameters for the shape-contrained smoothing function given passed to \code{smooth_method}. See help of the functions for the list of parameters.
#' @param alpha.control A list containing the controlling parameters for the
#'    alpha optimization steps of the algorithm. See \code{\link{alpha.setup}}.
#' @param algo.control A list containing controlling parameters for the
#'    whole algorithm. See \code{\link{algo.setup}}.
#' @param keep.trace Logical. If TRUE, the result contains a 'trace' element
#'    giving the intermediate values of alpha coefficients and ridge functions
#'    at each step of the algorithm. Useful for behaviour tracking.
#'
#' The CGAIM is expressed 
#'  \deqn{y_{i} = \beta_{0} + \sum_{j} \beta_{j} g_{j}(\alpha_{j}^{T} x_{ij})
#'    + \sum_{k} \gamma_{k} f_{k}(x_{ik}) + e_{i}}
#'  The formula interface considers \code{\link{g}} to identify index terms,
#'   \code{\link{s}} for smooth functions and can also include 
#'    linear terms. All smooth terms can be
#'    shape constrained.
#'
#' The CGAIM allows for linear constraints on the alpha coefficients. 
#'  Common constraints can be given for each index through the 
#'  \code{\link{g}} terms in the formula. Custom constraints, including multi-index constraints can be passed as a matrix through \code{alpha.control$Cmat}. See 
#'  \code{\link{alpha.setup}}.
#'
#' The CGAIM is fitted through an iterative algorithm that alternates between
#'  estimating the ridge functions \eqn{g_{j}} (and other non-index terms) and 
#'  updating the coefficients \eqn{\alpha_{j}}. The smoothing of ridge functions
#'  currently supports three methods: \code{\link[scam]{scam}} (the default), \code{\link[cgam]{cgam}} and \code{\link[scar]{scar}}. The list \code{smooth.control} allows controlling the functions.
#'
#' Updating the coefficient is made through quadratic programming. Currently,
#'  it is carried out by either the function 
#'  \code{\link[osqp:solve_osqp]{solve_osqp}} or 
#'  \code{\link[quadprog:solve.QP]{solve.QP}}. Controlling the updates is
#'  made through the parameter \code{alpha.control}. See \code{\link{alpha.setup}}
#'  for the detail.
#'
#' @return A \code{cgaim} object, i.e. a list with components:
#'  \item{alpha}{A named list of index coefficients.}
#'  \item{gfit}{A matrix containing the ridge and smooth functions 
#'    evaluated at the observations. Note that column ordering puts indices first and covariates after.}
#'  \item{indexfit}{A matrix containing the indices evaluated at the 
#'    observations.}
#'  \item{beta}{A vector containing the intercept and the scale coefficient
#'    of each ridge and smooth function. Includes the \eqn{\gamma_{k}} of
#'    the CGAIM model above. Note that ordering puts indices first and covariates after.}.
#'  \item{index}{A vector identifying to which index the columns of the
#'    element \code{x} belong.}
#'  \item{fitted}{A vector of fitted y values.}
#'  \item{residuals}{A vector of residuals.}
#'  \item{avcov}{The covariance matrix of alpha coefficients. Naively 
#'    computed from the last step of the algorithm. STILL EXPERIMENTAL.}
#'  \item{bvcov}{The covariance matrix of beta coefficients. 
#'    Computed from the design matrix created by ridge and smooth functions.
#'    STILL EXPERIMENTAL.}
#'  \item{gse}{Standard errors of ridge and smooth functions. Only
#'    available when \code{smooth_method = 'scam'}.}
#'  \item{rss}{The residual sum of squares of the fit.}
#'  \item{flag}{A flag indicating how the algorithm stopped. 1 for proper 
#'    convergence, 2 when the algorithm stopped for lack of step in the
#'    descent direction and 3 when the maximum number of iterations has
#'    been reached.}
#'  \item{trace}{If \code{keep.trace = TRUE}, a tracking list containing,
#'    for each iteration, values of alpha coefficients, 
#'    ridge and smooth function, convergence criteria, step length,
#'    as well as the number of iterations.}
#'  \item{covariates}{A data.frame containing the values of the variables not 
#'    entering any index.}
#'  \item{x}{A matrix containing the variables entering the indices. 
#'    The variables are mapped to each index throught the element \code{index}.}
#'  \item{y}{The response vector.}
#'  \item{weights}{The weights used for estimation.}
#'  \item{smooth.control}{The controlling parameters passed to the function
#'    for the smoothing steps.}
#'  \item{alpha.control}{The controlling parameters passed to the function
#'    for the alpha update steps.}
#'  \item{algo.control}{The controlling parameters passed to the function
#'    for the minimization algorithm.}
#'  \item{terms}{A \code{\link[stats]{terms.object}} representing the model.
#'    Useful for prediction.}
#'
#' @seealso \code{\link{confint.cgaim}} for confidence interval,
#'    \code{\link{predict.cgaim}} to predict new data,
#'    \code{\link{plot.cgaim}} to plot ridge function.
#'
#' @export
cgaim <- function(formula, data, weights, na.action, 
  smooth_method = c("scam", "cgam", "scar"), 
  smooth.control = list(), alpha.control = list(), algo.control = list(), 
  keep.trace = F)
{
  mt <- stats::terms(formula, specials = c("g", "s"), data = data)
  allvars <- all.vars(formula)
  gind <- attr(mt, "specials")$g
  p <- length(gind)
  ptot <- length(attr(mt, "term.labels"))
  if (missing(na.action)) na.action <- getOption("na.action")
  # Extract info on indices
  data <- stats::model.frame(stats::reformulate(allvars), data = data, 
    na.action = na.action)  
  index_interp <- stats::model.frame(mt[gind - 1], data = data)
  y <- stats::model.response(index_interp)
  attr(y, "varname") <- allvars[1]
  index_interp <- index_interp[-1]
  # Labels for smoothing step
  index_labels <- sapply(index_interp, attr, "label")
  ulabs <- unique(index_labels)
  plabs <- length(ulabs)
  if (plabs < p){
    for (ilab in ulabs){
      wlabs <- which(index_labels == ilab)
      if (length(wlabs) > 1){
        index_labels[wlabs] <- sprintf("%s%i", ilab, 1:length(wlabs))
      }
    }
  }
  # Organizing variables in indices as a matrix with index positions   
  pvec <- sapply(index_interp, attr, "nterms")
  index <- rep(1:p, pvec)
  names(index) <- rep(index_labels, pvec)
  n <- nrow(index_interp)
  Xind <- matrix(0, n, sum(pvec))
  for (j in 1:p){
    Xind[,index == j] <- index_interp[[j]]
    colnames(Xind)[index == j] <- colnames(index_interp[[j]])
  }
  # Initialize alpha controls
  alpha.control <- do.call(alpha.setup, alpha.control)
  nai <- length(alpha.control$alpha.start)
  if (nai > 0 && nai != sum(pvec)){
    warning("alpha.init length is inconsistent with index matrix and is recycled")
    alpha.control$alpha.start <- rep_len(alpha.control$alpha.start, sum(pvec))
  }
  # Organizing the constraint matrix
  gconsts <- as.matrix(Matrix::bdiag(lapply(index_interp, attr, "Cmat")))
  constr_mat <- alpha.control$Cmat
  if (!is.null(constr_mat)){
    if (ncol(constr_mat) != ncol(Xind)){
      stop("Number of columns in alpha.control$Cmat does not match
        the indices")
    }
  }
  if (length(gconsts) > 0 || length(constr_mat) > 0){
    Cmat <- rbind(constr_mat, gconsts)
    Cmat <- Cmat[!duplicated(Cmat),, drop = F]
    alpha.control$Cmat <- Cmat
  }
  # Prepare parameters for smoothing
  smooth_method <- match.arg(smooth_method)
  setup_fun <- sprintf("%s.setup", smooth_method)
  smooth.control <- do.call(setup_fun, list(mt, data, smooth.control))
  algo.control$smooth_method <- smooth_method
  # Calling the fitting function
  pars <- algo.control[names(algo.control) %in% 
    methods::formalArgs(gaim_gn)]
  pars$x <- Xind
  pars$y <- y
  if (missing(weights)){
    pars$w <- rep(1, n)
  } else {
    pars$w <- rep_len(weights, n)
  }
  pars$index <- index
  pars$smooth.control <- smooth.control
  pars$alpha.control <- alpha.control
  # Fitting the model
  result <- do.call(gaim_gn, pars)
  result$alpha <- split(result$alpha, index)
  for (j in 1:p){
    names(result$alpha[[j]]) <- colnames(index_interp[[j]])
  }
  names(result$alpha) <- index_labels
  attributes(result$gfit) <- attributes(result$gfit)[c("dim", "dimnames")]
  if (!is.null(result$bvcov)){
    colnames(result$bvcov) <- rownames(result$bvcov) <- names(result$beta)
  }  
  result$covariates <- smooth.control$Xcov
  result$x <- Xind
  result$y <- y
  result$weights <- pars$w
  result$smooth.control <- smooth.control
  result$alpha.control <- alpha.control
  result$algo.control <- algo.control[names(algo.control) %in% 
    methods::formalArgs(gaim_gn)]
  result$terms <- mt
  class(result) <- "cgaim"
  return(result)
}

gaim_gn <- function(x, y, w, index, 
  smooth.control = list(), alpha.control = list(),
  keep.trace = FALSE, max.iter = 50, tol = 1e-3, min.step.len = 0.1, 
  halving = T, convergence_criterion = c("change", "offset"), 
  smooth_method = "scam")
{
  convergence_criterion <- match.arg(convergence_criterion)
  if (convergence_criterion == "change") tol <- rep_len(tol, 2)
  # Useful objects
  n <- length(y)
  p <- max(index)
  d <- length(index)
  ind_pos <- split(1:d, index)
  # Alpha initialization
  if (is.null(alpha.control$alpha.start)){    
    init_pars <- alpha.control
    init_pars <- within(init_pars,{
      y = y; x = x; w = w; index = index
    })
    alpha <- do.call(alpha_init, init_pars)
  } else {
    alpha <- alpha.control$alpha.start
  }
  # Indice computation
  zs <- sapply(ind_pos, function(i, x, a) x[,i] %*% a[i], 
    x = x, a = alpha)
  colnames(zs) <- unique(names(index))
  # Initial smoothing
  smooth_fun <- sprintf("smooth_%s", smooth_method) 
  smo_par <- c(smooth.control, 
    list(y = y, x = zs, weights = w))
  gz <- do.call(smooth_fun, smo_par)   
  yhat <- gz$intercept + rowSums(gz$gz)
  r <- y - yhat  
  # Convergence criterion
  eps <- tol + 1
  c1 <- 1      
  l2 <- L2(y, yhat, w)  
  # If requested: tracing the algorithm evolution
  if (keep.trace){
    trace.list <- list(
      criterion = matrix(NA, max.iter, length(tol)),
      alpha = matrix(NA, nrow = max.iter, ncol = d),
      gfit = rep(list(NA), max.iter),
      step.len = rep(NA, max.iter) 
    )
    trace.list$criterion[1,] <- eps
    trace.list$step.len[1] <- 1  
    trace.list$alpha[1,] <- alpha 
    trace.list$gfit[[1]] <- gz$gz
  }
  # Gauss-Newton search
  alpha_pars <- c(
    alpha.control[names(alpha.control) %in% 
      methods::formalArgs(alpha_update)],
    list(x = x, w = w, index = index, delta = TRUE))
  stopflag <- 0
  while(stopflag == 0){
    # Update alphas
    alpha_pars$y <- r
    alpha_pars$dgz <- gz$dgz
    alpha_pars$alpha <- alpha
    delta <- do.call(alpha_update, alpha_pars)
    # Halving in case of bad steps
    cuts <- 1
    repeat {   
      delta <- delta * cuts
      alpha.new <- alpha + delta
      alpha.new <- unlist(
        tapply(alpha.new, index, normalize, type = alpha.control$norm.type))
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
    if (convergence_criterion == "change"){
      eps <- c((l2 - l2.new) / l2, max(abs(alpha.new - alpha) / abs(alpha)))
    } else {
      eps <- offset_convergence(r, x, gz$dgz[,index])
    }
    r <- y - yhat
    alpha <- alpha.new
    l2 <- l2.new
    c1 <- c1 + 1       
    if (keep.trace){ # Tracing the algorithm evolution
      trace.list$criterion[c1,] <- eps
      trace.list$step.len[c1] <- cuts 
      trace.list$alpha[c1,] <- alpha 
      trace.list$gfit[[c1]] <- gz$gz
    }
    if (all(eps < tol)){
      stopflag <- 1
    } else { 
      if (c1 >= max.iter){
        stopflag <- 3
        warning(sprintf("Fitting did not converge after %i iterations. Consider revising the constraints", max.iter))
      } 
    }    
  }
  final.gz <- scale(gz$gz)
  final.gz[,apply(gz$gz, 2, stats::sd) == 0] <- 
    gz$gz[,apply(gz$gz, 2, stats::sd) == 0]
  betas <- attr(final.gz, "scaled:scale")
  r <- y - yhat
  sig2 <- sum(r^2) / (n - d)
  # Parameter covar matrix
  Vmat <- x * gz$dgz[,index]
  vtv <- crossprod(Vmat)
  avcov <- try(chol2inv(chol(vtv)), silent = TRUE)
  if (!inherits(avcov, "try-error")){
    avcov <- avcov * sig2
    colnames(avcov) <- rownames(avcov) <- 
      paste(names(index), colnames(x), sep = ".")
  } else {
    avcov <- NULL
  }
  bvtv <- crossprod(cbind(1, final.gz))
  bvcov <- try(chol2inv(chol(bvtv)), silent = TRUE)
  if (inherits(bvcov, "try-error")){
    bvcov <- NULL
  }
  output <- list(alpha = alpha, gfit = final.gz, indexfit = zs, 
    beta = c(gz$intercept + sum(attr(final.gz, "scaled:center")), betas),
    index = index, fitted = yhat, residuals = r, avcov = avcov, bvcov = bvcov,
    gse = gz$se, rss = l2, flag = stopflag, niter = c1)
  if (keep.trace){
    if (convergence_criterion == "change"){
      colnames(trace.list$criterion) <- c("rss", "alpha")
    } else {
      colnames(trace.list$criterion) <- "offset"
    }
    trace.list$criterion <- trace.list$criterion[1:c1,,drop = FALSE]
    trace.list$alpha <- trace.list$alpha[1:c1,]
    trace.list$gfit <- trace.list$gfit[1:c1]
    trace.list$step.len <- trace.list$step.len[1:c1]
    output$trace <- trace.list
  }
  return(output)   
}