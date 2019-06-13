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
  if (any(sapply(x, nrow) != n)) stop("Matrices in x must have number of rows equal to the length of y")
  p <- length(x)
  pvec <- sapply(x, ncol)
  if (is.null(names(x))) names(x) <- sprintf("V%i", 1:p)
  xnames <- names(x)
  if(missing(w)) w <- rep(1/n,n) # Weights
  # Default values for controlling the algorithm
  defalgo.control <- list(type = "gauss.newton")
  algo.control <- c(algo.control, 
    defalgo.control[!names(defalgo.control) %in% names(algo.control)])
  # Default values for controlling alpha updates
  defalpha.control <- list(norm.type = "L2", monotone = 0, sign.const = 0, 
    init.type = "regression", constraint.algo = "QP", 
    delta = FALSE)
  alpha.control <- c(alpha.control, 
    defalpha.control[!names(defalpha.control) %in% names(alpha.control)])
  alpha.control[c("norm.type", "sign.const", "monotone")] <- 
    lapply(alpha.control[c("norm.type", "sign.const", "monotone")], 
    rep_len, length.out = p)
  defsmoo.control <- list(shape = "tp")
  smooth.control <- c(smooth.control, 
    defsmoo.control[!names(defsmoo.control) %in% names(smooth.control)])
  smooth.control$shape <- rep_len(smooth.control$shape, length.out = p)
  # Center y
  beta0 <- mean(y)
  y <- y - beta0
  # Initialize alphas
  gn_par <- c(alpha.control[names(alpha.control) %in% names(formals(gn_update))], #! Pas nécessaire?
    list(r = y, x = x, dgz = rep(list(1), p), w = w))
  alpha <- do.call(gn_update, gn_par)
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
  result$fitted <- yhat
  return(result)
}

#' AIM through genetic algorithm
aim.ga <- function(x, y, w, alpha = rep(0, sum(pvec)), trace = FALSE, 
   smooth.control = list(), alpha.control = list(), algo.control = list()) 
{
  p <- length(x)
  pvec <- sapply(x, ncol)
  pind <- rep(1:p, pvec)
  alphavec <- unlist(alpha)
  f <- function(alphavec, x, y, w){
    alpha <- split(alphavec, pind)
    alpha <- Map(normalize, alpha, alpha.control$norm.type)
    zs <- mapply("%*%", x, alpha)
    smo_par <- c(smooth.control, list(y = y, x = zs, w = w))
    gz <- do.call(smoothing, smo_par)
    yhat <- gz$intercept + rowSums(gz$gz)
    l2 <- L2(y, yhat, w)
    return(l2)
  }
  monotone.cons <- function(alphavec){
    consts <- const.matrix(pind, alpha.control$monotone, rep(0, p))
    consts %*% alphavec
  }
  objfun <- function(alphavec, x, y, w){
    pen <- sqrt(.Machine$double.xmax)
    constvec <- monotone.cons(alphavec)
    penalty <- apply(-constvec, 1, max, 0) * pen
    -f(alphavec, x, y, w) - sum(penalty)
  }
  ga_par <- algo.control[names(algo.control) %in% names(formals(ga))]
  if (is.null(alpha.control$upper)){
    alpha.control$upper <- lapply(pvec, rep, x = 1)
  }
  if (is.null(alpha.control$lower)){
    alpha.control$lower <- lapply(pvec, rep, x = -1)
  }
  lower <- unlist(alpha.control$lower)
  lower[alpha.control$sign.const[pind] == 1] <- 0
  upper <- unlist(alpha.control$upper)
  upper[alpha.control$sign.const[pind] == -1] <- 0
  ga_par <- within(ga_par, {
    type <- "real-valued" 
    fitness <- objfun
    x <- x
    y <- y
    w <- w
    suggestions <- alphavec
    lower <- lower
    upper <- upper
    monitor <- trace
  })
  result <- do.call(ga, ga_par)  
  new.alphavec <- result@solution[1,]
  newalpha <- split(new.alphavec, pind)  
  zs <- mapply("%*%", x, newalpha)
  smo_par <- c(smooth.control, list(y = y, x = zs, w = w))
  gz <- do.call(smoothing, smo_par)
  output <- list(alpha = newalpha, gz = gz$gz, z = zs) 
}

#' 
aim.optim <- function(x, y, w, alpha = rep(0, sum(pvec)), 
   smooth.control = list(), alpha.control = list(), algo.control = list()) 
{
  p <- length(x)
  pvec <- sapply(x, ncol)
  pind <- rep(1:p, pvec)
  alphavec <- unlist(alpha)
  f <- function(alphavec, x, y, w){
    alpha <- split(alphavec, pind)
    alpha <- Map(normalize, alpha, alpha.control$norm.type)
    zs <- mapply("%*%", x, alpha)
    smo_par <- c(smooth.control, list(y = y, x = zs, w = w))
    gz <- do.call(smoothing, smo_par)
    yhat <- gz$intercept + rowSums(gz$gz)
    l2 <- L2(y, yhat, w)
    return(l2)
  }
  constr <- with(alpha.control, all(monotone != 0) & all(sign.const != 0))
  if (constr){
    ui <- const.matrix(pind, alpha.control$monotone, alpha.control$sign.const)
    ui <- t(ui)
    ci <- rep(0, nrow(ui))
    adjust <- ui %*% alphavec - ci <= 0
    while(any(adjust)){
      alphavec <- alphavec + 
        t(ui[adjust,,drop = F]) %*% rep(.001, sum(adjust))
      adjust <- ui %*% alphavec - ci <= 0
    }    
    optim_par <- c(
      algo.control[names(algo.control) %in% names(formals(constrOptim))],
      list(theta = alphavec, f = f, grad = NULL, ui = ui, ci = ci, 
        x = x, y = y, w = w))
    result <- do.call(constrOptim, optim_par)
  } else {
    optim_par <- c(
      algo.control[names(algo.control) %in% names(formals(constrOptim))],
      list(par = alphavec, fn = f, x = x, y = y, w = w))
    result <- do.call(optim, optim_par)
  }
  new.alphavec <- result$par
  newalpha <- split(new.alphavec, pind)  
  zs <- mapply("%*%", x, newalpha)
  smo_par <- c(smooth.control, list(y = y, x = zs, w = w))
  gz <- do.call(smoothing, smo_par)
  output <- list(alpha = newalpha, gz = gz$gz, z = zs) 
}

#' AIM through backfitting of single-index models
aim.backfit <- function(x, y, w, alpha = rep(0, sum(pvec)), gz, 
   smooth.control = list(), alpha.control = list(), algo.control = list(), 
   trace = T) 
{
    n <- length(y)
    p <- length(x)
    pvec <- sapply(x, ncol)
    def.control <- list(bf.tol = 5e-03, bf.maxit = 50)
    algo.control <- c(algo.control, 
      def.control[!names(def.control) %in% names(algo.control)]) 
    if (missing(gz)){
      gz <- list(intercept = mean(y), gz = matrix(0, n, p), 
        dgz = matrix(1, n, p))
    }
    gs <- gz$gz
    if (missing(alpha)) alpha <- sapply(pvec, rep, x = 0)
    # Backfitting
    rel.delta <- algo.control$bf.tol + 1
    cbf <- 0
    if (trace){ # Tracing the algorithm evolution
        trace.list <- list(
            criteria = matrix(NA, nrow = algo.control$bf.maxit, ncol = 3, 
              dimnames = list(NULL, c("delta", "rel.delta", "rss"))), 
            alpha = mapply(matrix, ncol = pvec, 
              MoreArgs = list(nrow = algo.control$bf.maxit)), 
            gz = replicate(p, matrix(nrow = n, ncol = algo.control$bf.maxit), 
              simplify = F)
        )
        names(trace.list$gz) <- names(x)
    }
    while (rel.delta > algo.control$bf.tol && cbf < algo.control$bf.maxit){
        deltaf <- 0
        f0norm <- sum(apply(gs^2, 2, weighted.mean, w = w))
        for (j in 1:p){ #! /!\ Oscillation of betas and gz
            rj <- y - rowSums(gs[,-j])
            al.co.j <- alpha.control
            al.co.j[c("norm.type", "monotone", "sign.const")] <- 
              sapply(al.co.j[c("norm.type", "monotone", "sign.const")], "[", j)
            sm.co.j <- smooth.control
            sm.co.j$type <- sm.co.j$type[j]
            gzj <- within(gz,{
              gz <- gz[,j, drop = F]
              dgz <- dgz[,j, drop = F]
            })
            jterm <- aim.GaussNewton(x = x[j], y = rj, w = w, 
              alpha = alpha[j], gz = gzj, 
              smooth.control = sm.co.j, alpha.control = al.co.j,
              algo.control = algo.control, trace = trace)
            alpha[[j]] <- jterm$alpha[[1]]
            deltaf <- deltaf + weighted.mean((jterm$gz - gs[,j])^2, w)  
            gs[,j] <- jterm$gz
        }
        rel.delta <- sqrt(deltaf/f0norm) 
        cbf <- cbf + 1
        if (trace){ # Tracing the algorithm evolution
          trace.list$criteria[cbf,] <- c(deltaf, rel.delta, 
            sum((y - rowSums(gs[,-j]))^2))
          for (i in 1:p){
            trace.list$alpha[[i]][cbf,] <- alpha[[i]]
            trace.list$gz[[i]][,cbf] <- gs[,i]
          }              
        }
    }
    if (cbf == algo.control$bf.maxit) warning(sprintf("Convergence not attained after %i iterations", algo.control$bf.maxit))
    zs <- mapply("%*%", x, alpha)
    output <- list(alpha = alpha, gz = gs, z = zs)
    if (trace) output$trace <- trace.list
    return(output)
}

#' AIM by Gauss-Newton algorithm
#'
#' Repeat two steps until convergence:
#'    1) Update coefficient of all indices at once by Gauss-Newton
#'    2) Apply GAM (or SCAM) on indices to update Ridge Functions.
aim.GaussNewton <- function(x, y, w, alpha, gz, 
  smooth.control = list(), alpha.control = list(), algo.control = list(), 
  trace = FALSE)
{
    n <- length(y)
    p <- length(x)
    pvec <- sapply(x, ncol)
    d <- sum(pvec)
    pind <- rep(1:p, pvec)
    defalgo.control <- list(tol = 0.001, max.iter = 50, min.step.len = 0.1, 
      halving = T, convergence_criterion = "LS")
    algo.control <- c(algo.control,defalgo.control[!names(defalgo.control) %in% 
      names(algo.control)])
#    defalpha.control <- list(init.type = "runif", init.first.pos = TRUE)
#    alpha.control <- c(alpha.control, 
#      defalpha.control[!names(defalpha.control) %in% names(alpha.control)])
    if (missing(gz)){
      gz <- list(intercept = mean(y), gz = matrix(0, n, p), 
        dgz = matrix(1, n, p))
    }
    if (missing(alpha)) alpha <- sapply(pvec, rep, x = 0)
    # Convergence criterion    
    yhat <- gz$intercept + rowSums(gz$gz)
    r <- y - yhat
    l2 <- L2(y, yhat, w)
    eps <- algo.control$tol + 1
     # Tracing the algorithm evolution
    if (trace){
      trace.list <- list(
        convergence = rep(NA, algo.control$max.iter),
        alpha = lapply(pvec, matrix, data = NA, nrow = algo.control$max.iter),
        gz = array(NA, dim = c(n, p, algo.control$max.iter))
      )
      trace.list$L2[1] <- l2 
      for (j in 1:p) trace.list$alpha[[j]][1,] <- alpha[[j]] 
      trace.list$gz[,,1] <- gz$gz
    }
    # Gauss-Newton search
    c1 <- 1
    gn_par <- c(
      alpha.control[names(alpha.control) %in% names(formals(gn_update))],
      list(x = x, w = w))
    smo_par <- c(smooth.control, list(y = y, w = w))
    while(eps > algo.control$tol && c1 <= algo.control$max.iter){
        gn_par$r <- r
        gn_par$dgz <- gz$dgz
        gn_par$alpha <- alpha
        # Update alphas
        #! Work with alphas as vectors (less operations)
        delta <- do.call(gn_update, gn_par)
        cuts <- 1
        l2.old <- l2
        repeat {   
        #Loop for halving steps in case of bad step (backtracking line search)
            delta <- lapply(delta, "*", cuts)
            alpha.new <- Map("+", alpha, delta)
            z <- mapply("%*%", x, alpha.new)
            smo_par$x <- z
            gz <- do.call(smoothing, smo_par)
            yhat <- gz$intercept + rowSums(gz$gz)
            l2 <- L2(y, yhat, w)
            if (l2 < l2.old || cuts < algo.control$min.step.len || 
              !algo.control$halving){
               break
            }
            cuts <- cuts / 2
        }
        r <- y - yhat
        alpha.new <- Map(normalize, alpha.new, alpha.control$norm.type)
        eps <- switch(algo.control$convergence_criterion,
          LS = (l2.old - l2)/l2.old,
          alpha = max(abs(unlist(alpha.new) - unlist(alpha)) / 
            abs(unlist(alpha))),
          orthogonal = { #! Doesn't work well
            Q <- qr.qty(qr(
              Reduce("cbind", Map("*", x, as.data.frame(gz$dgz)))), r)
            (sqrt(sum(Q[1:d]^2)) / sqrt(d)) / 
              (sqrt(sum(Q[-(1:d)]^2)) / sqrt(n - d))
          }
        )
        alpha <- alpha.new       
        if (trace){ # Tracing the algorithm evolution
          trace.list$convergence[c1] <- l2 
          for (j in 1:p) trace.list$alpha[[j]][c1,] <- alpha[[j]] 
          trace.list$gz[,,c1] <- gz$gz
        }
        c1 <- c1 + 1
    }
    z <- mapply("%*%", x, alpha)
    smo_par$x <- z
    gz <- do.call(smoothing, smo_par)
    output <- list(alpha = alpha, gz = gz$gz, z = z)
    if (trace) output$trace <- trace.list
    return(output)   
}