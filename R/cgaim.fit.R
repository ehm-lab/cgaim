#' @export
cgaim.fit <- function(x, y, w, index, 
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
    niter = c1, edf = edf, gcv = gcv, dg = gz$dgz, gse = gz$se, 
    active = alphaup$active)
  if (trace){
    trace.list$criterion <- trace.list$criterion[1:c1]
    trace.list$alpha <- trace.list$alpha[1:c1,]
    trace.list$gfit <- trace.list$gfit[1:c1]
    trace.list$step.len <- trace.list$step.len[1:c1]
    output$trace <- trace.list
  }
  return(output)   
}