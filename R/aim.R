#!############################################################################
#!
#!           A global AIM
#!
#!############################################################################

#' Additive index models
#'
#' Fit an additive index model (AIM). Possibility to use different algorithms.
#'    Still a test function.
#'
#' @param x A list of matrices giving the data for each index. If a matrix is
#'    given, a single-index model is fitted.
#' @param y A numeric vector containing the output of the model.
#' @param w A numeric vector containing weights.
#' @param smooth.control A list containing the controlling parameters for the
#'    smoothing of the ridge functions g.
#' @param alpha.control A list containing the controlling parameters for the
#'    estimation of the index coefficients alpha.
#' @param algo.control A list containing the controlling parameters for the
#'    algorithm. Notably includes the type of algorithm to fit the AIM.
#' @param trace Logical indicating if the algorithm should be traced.
aim <- function(y, x, w, smooth.control = list(), alpha.control = list(),  
  algo.control = list(), trace = FALSE)
{
  y <- as.vector(y)
  n <- length(y)
  if (!is.list(x)) x <- list(x) #! data.frame
  x <- lapply(x, as.matrix)
  if (any(sapply(x, nrow) != n)) stop("Matrices in x must have number of rows equal to the length of y")
  p <- length(x)
  pvec <- sapply(x, ncol)
  indices <- pvec > 1
  pdex <- sum(indices)
  if (is.null(names(x))) names(x) <- sprintf("V%i", 1:p)
  xnames <- names(x)
  if(missing(w)) w <- rep(1/n,n) # Weights
  # Default values for controlling the algorithm
  defalgo.control <- list(tol = 5e-3, max.iter = 50, min.step.len = 0.1, 
      halving = T, convergence_criterion = "LS")
  algo.control <- c(algo.control, 
    defalgo.control[!names(defalgo.control) %in% names(algo.control)])
  # Default values for controlling alpha updates
  defalpha.control <- list(norm.type = "L2", monotone = 0, sign.const = 0, 
    init.type = "regression", constraint.algo = "QP", 
    delta = TRUE)
  alpha.control <- c(alpha.control, 
    defalpha.control[!names(defalpha.control) %in% names(alpha.control)])
  alpha.control[c("norm.type", "sign.const", "monotone")] <- 
    lapply(alpha.control[c("norm.type", "sign.const", "monotone")], 
    rep_len, length.out = pdex)
  defsmoo.control <- list(shape = "tp")
  smooth.control <- c(smooth.control, 
    defsmoo.control[!names(defsmoo.control) %in% names(smooth.control)])
  smooth.control$shape <- rep_len(smooth.control$shape, length.out = p)
  # Center y
  beta0 <- mean(y)
  y <- y - beta0
  # Initialize alphas
  gn_par <- c(alpha.control[names(alpha.control) %in% names(formals(gn_update))], #! Pas nécessaire?
    list(r = y, x = x[indices], dgz = rep(list(1), pdex), w = w))
  gn_par$delta <- FALSE
  alpha <- rep(list(1), p)
  alpha[indices] <- do.call(gn_update, gn_par)
  # Initialize ridge functions gz
  zs <- mapply("%*%", x, alpha)
  smo_par <- c(smooth.control, list(y = y, x = zs, w = w))
  gz <- do.call(smoothing, smo_par)
  # Main algorithm
  result <- switch(algo.control$type,
    two.steps = list(alpha = alpha, gz = gz$gz, z = zs),
    gauss.newton = aim.GaussNewton(x = x, y = y, w = w, alpha = alpha, gz = gz,
      smooth.control = smooth.control, alpha.control = alpha.control, 
      algo.control = algo.control, trace = trace),
    backfitting = aim.backfit(x = x, y = y, w = w, alpha = alpha, gz = gz, 
      smooth.control = smooth.control, alpha.control = alpha.control, 
      algo.control = algo.control, trace = trace),
    optim = aim.optim(x = x, y = y, w = w, alpha = alpha, 
      smooth.control = smooth.control, alpha.control = alpha.control, 
      algo.control = algo.control),
    ga = aim.ga(x = x, y = y, w = w, alpha = alpha, 
      smooth.control = smooth.control, alpha.control = alpha.control, 
      algo.control = algo.control, trace = trace),
    stop("Unknown algo type")
  )  
  names(result$alpha) <- names(x)
  colnames(result$gz) <- names(x)
  # Adjust final values
  final.gz <- scale(result$gz)
  betas <- attr(final.gz, "scaled:scale")
  yhat <-  beta0 + final.gz %*% betas
  result$gz <- final.gz
  result$coef <- c(beta0, betas)
  names(result$coef) <- c("intercept", names(x))
  result$fitted <- drop(yhat)
  return(result)
}



#' AIM by Gauss-Newton algorithm
#'
#' x is a matrix with the index variables
#' y is the response
#' w are weigths
#' index is a vector identifying to which index each variable in x belongs. 
#'      Additional variables not part of any index are not represented in index.
#'      Overall we consider that the first d variables are in indices and
#'      the d + 1 to dtot are additional covariates
#'
#' Repeat two steps until convergence:
#'    1) Update coefficient of all indices at once by Gauss-Newton
#'    2) Apply GAM (or SCAM) on indices to update Ridge Functions.
gaim_gn <- function(x, y, w, index, 
  smooth.control = list(), alpha.control = list(),
  keep.trace = FALSE, max.iter = 50, tol = 1e-3, min.step.len = 0.1, 
  halving = T, convergence_criterion = c("LS", "Alpha", "Offset"))
{
  convergence_criterion <- match.arg(convergence_criterion)
  # Useful objects
  n <- length(y)
  p <- max(index)
  d <- length(index)
  dtot <- ncol(x) 
  ind_pos <- split(1:d, index)
  xind <- x[,1:d]
  if (dtot > d){
    xadd <- x[,(d + 1):dtot]
  } else {
    xadd <- NULL
  }     
  # Alpha initialization    
  init_pars <- alpha.control[names(alpha.control) %in%
    intersect(names(formals(alpha_init)), names(formals(alpha_update)))
  ]
  init_pars <- within(init_pars,{
    y = y; x = xind; w = w; index = index
  })
  alpha <- do.call(alpha_init, init_pars)
  # Indice computation
  zs <- sapply(ind_pos, function(i, x, a) x[,i] %*% a[i], 
    x = xind, a = alpha)
  xsmoo <- cbind(zs, xadd)
  # Initial smoothing 
  smo_par <- c(smooth.control, list(y = y, x = xsmoo, w = w))
  gz <- do.call(smoothing, smo_par)   
  yhat <- gz$intercept + rowSums(gz$gz)
  r <- y - yhat  
  # Convergence criterion
  eps <- tol + 1
  c1 <- 1      
  l2 <- L2(y, yhat, w)  
  # If requested: tracing the algorithm evolution
  if (keep.trace){
    trace.list <- list(
      convergence = rep(NA, max.iter),
      alpha = matrix(NA, nrow = max.iter, ncol = d),
      gz = rep(list(NA), max.iter)
    )
    trace.list$convergence[1] <- eps 
    trace.list$alpha[1,] <- alpha 
    trace.list$gz[[1]] <- gz$gz
  }
  # Gauss-Newton search
  alpha_pars <- c(alpha.control[names(alpha.control) %in%
      names(formals(alpha_update))],
    list(x = xind, w = w, index = index))
  alpha_pars$delta <- TRUE
  while(eps > tol && c1 < max.iter){
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
        x = xind, a = alpha.new)
      smo_par$x <- cbind(zs, xadd)
      gz <- do.call(smoothing, smo_par)
      yhat <- gz$intercept + rowSums(gz$gz)
      l2.new <- L2(y, yhat, w)
      r <- y - yhat
      eps <- switch(convergence_criterion,
        LS = (l2 - l2.new) / l2,
        Alpha = max(abs(alpha.new - alpha) / abs(alpha)),
        Offset = offset_convergence(r, xind, gz$dgz[,index])
      )
      if (eps > 0 || cuts < min.step.len || 
        !halving){
         break
      }
      cuts <- cuts / 2
    }
    alpha <- alpha.new
    l2 <- l2.new
    c1 <- c1 + 1       
    if (keep.trace){ # Tracing the algorithm evolution
      trace.list$convergence[c1] <- eps 
      trace.list$alpha[c1,] <- alpha 
      trace.list$gz[[c1]] <- gz$gz
    }
  }
  if (keep.trace){
    trace.list$convergence <- trace.list$convergence[1:c1]
    trace.list$alpha <- trace.list$alpha[1:c1,]
    trace.list$gz <- trace.list$gz[1:c1]
    trace.list$niterations <- c1
  }
  final.gz <- scale(gz$gz)
  final.gz[,apply(gz$gz, 2, sd) == 0] <- gz$gz[,apply(gz$gz, 2, sd) == 0]
  betas <- attr(final.gz, "scaled:scale")
  output <- list(alpha = alpha, gz = final.gz, z = zs, 
    beta = c(gz$intercept, betas), am.fit = gz$fit,
    x = x, y = y, index = index, dgz = gz$dgz, fitted = yhat)
  if (keep.trace) output$trace <- trace.list
  return(output)   
}