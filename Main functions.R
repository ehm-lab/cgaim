#!#############################################################################
#!
#!                        Main GAIM functions                                  
#!
#!#############################################################################

source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/Package functions/Secondary functions.R")

#! fit single term ppr
single.term <- function(x, y, w, alpha = alpha.init(ncol(x)), beta1 = 1, smoother = "spline", smoother.args = NULL, control = NULL) 
# y: vector of response values
# x: matrix of predictor values
# w: weights
# beta1: initial value of beta1 (possibly inherited from previous passes by the backfitting)
# smoother: the smoother used
# control: control parameters for the search algorithm
{
    def.control <- list(tol = 5e-3, max.iter = 20, min.step.len = 0.1, halving = T)
    control <- c(control,def.control[!names(def.control) %in% names(control)])
    smooth.fun <- sprintf("%s.smoother", smoother)
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    if(missing(w)) w <- rep(1/n,n)
    # First iteration
#    alpha <- alpha.init(type = "ols", p = p, y = y, x = x, normalize = TRUE)
    z <- x %*% alpha 
    g <- do.call(smooth.fun, c(list(y = y/beta1, z = z, w = w), smoother.args))
#    x11(); plot(x%*%alpha, g$gz)
#    z <- x %*% alpha
#    sm.fit <- do.call("smooth.spline", list(x = z, y = y/beta1))
#    gz <- predict(sm.fit, z)$y
#    gznorm <- sqrt(sum(g$gz^2))
#    g$gz <- g$gz / gznorm
#    beta1 <- beta1*gznorm
    l2 <- L2(y, as.vector(beta1*g$gz), w)
    eps <- (var(y) - l2)/var(y)
    # Iterate until convergence
    c1 <- 1
    while (eps > control$tol && c1 <= control$max.iter){        
#        dgz <- do.call(deriv.fun, list(x = x, z = z, gz = gz, alpha = alpha)) #! TODO: Modify to allow passing function for derivation<
#        dgz <- predict(sm.fit, z, deriv = 1)$y 
        # As in Roosen and Hastie (1994)
        r <- as.vector(y) - g$gz
        xf <- beta1 * x * matrix(g$dgz, nrow = n, ncol = p, byrow = F)
        delta <- coef(lm(r ~ 0 + xf, weights = w))
#        delta <- apply(beta1*x, 2, function(x){
#            xf <- x*g$dgz
#            coef(lm(r ~ 0 + xf, weights = w))
#        })
        cuts <- 1
        l2.old <- l2
        repeat {   #Loop for halving steps in case of bad step  #
            delta <- delta * cuts
            alpha.new <- alpha + delta
#            alpha.new <- sign(alpha.new[1])*alpha.new/sqrt(sum(alpha.new^2)) # For identifiability
            z <- x %*% alpha.new
            g <- do.call(smooth.fun, c(list(y = y/beta1, z = z, w = w), smoother.args))
#            z <- x %*% alpha.new
#            sm.fit <- do.call(smoother, list(x = z, y = y/beta1))
#            gz <- predict(sm.fit, z)$y
#            gznorm <- sqrt(sum(g$gz^2))
#            g$gz <- g$gz / gznorm
#            beta1 <- beta1*gznorm
            l2 <- L2(as.vector(y), beta1*g$gz, w)
            if (l2 < l2.old || cuts < control$min.step.len || !control$halving){
               break
            }
            cuts <- cuts / 2
        }
        alpha <- alpha.new
        eps <- (l2.old - l2)/l2.old
        c1 <- c1 + 1
#        plot(x%*%alpha, g$gz)
#        print(alpha)
    }
    alpha <- alpha/sqrt(sum(alpha^2)) 
    z <- x %*% alpha
    g <- do.call(smooth.fun, c(list(y = y/beta1, z = z, w = w), smoother.args))
#    sm.fit <- do.call(smoother, list(x = z, y = y/beta1))
#    gz <- predict(sm.fit, z)$y
    gznorm <- sqrt(sum(g$gz^2))
    g$gz <- g$gz / gznorm
    beta1 <- beta1*gznorm
    return(list(alpha = alpha, z = z, gz = g$gz, beta = beta1)) #! TODO: add predict
}

#! Additive index models
aim <- function(x, y, w, smoother = "spline", smoother.args = list(NULL), control = NULL, init.type = "runif", trace = T) 
# x: a list of matrices giving the variables for each index. The matrices need have the same number of rows (individuals) but can have diffeent numbers of columns (variables). If a matrix is provided, a single-index model is fitted. 
# y: a numeric vector containing the response observations
# w: weigths
# smoother: a character vector giving the smoother used for each index. Can be different for each index. Recycled if necessary.
# smoother.args: a (nested) list containing optional arguments for each smoother. 
# control: control parameters for the search algorithm
# trace: if true, keep criteria and parameters values of each iteration
{
    #! Estimate all coefficients first? As in Wang et al. (2015)
    if (!is.list(x)) x <- list(x)
    xnames <- names(x)
    n <- length(y)
    p <- length(x)
    pvec <- sapply(x, ncol)
    if (any(sapply(x, nrow) != n)) stop("Matrices in x must have number of rows equal to the length of y")
    smoother <- rep_len(smoother, p)
    if (length(smoother.args) != p){
       warning("Smoother.args does not have the same length as x and is thus recycled.")
       smoother.args <- rep_len(smoother.args, p)
    }
    if(missing(w)) w <- rep(1/n,n)
    def.control <- list(bf.tol = 5e-03, bf.maxit = 50)
    control <- c(control,def.control[!names(def.control) %in% names(control)]) 
    # Initialize variables
    beta0 <- mean(y) #! Possibly remove for the fisher scoring
    betas <- rep(1,p)
    alphas <- mapply(alpha.init, p = pvec, x = x, MoreArgs = list(y = y, type = init.type))
    gs <- matrix(0, n, p)
    # Backfitting
    rel.delta <- control$bf.tol + 1
    cbf <- 0
    if (trace){ # Tracing the algorithm evolution
        trace.list <- list(
            criteria = matrix(NA, nrow = control$bf.maxit, ncol = 3, dimnames = list(NULL, c("delta", "rel.delta", "rss"))), 
            beta = matrix(NA, nrow = control$bf.maxit, ncol = p, dimnames = list(NULL, xnames)), 
            alpha = mapply(matrix, ncol = pvec, MoreArgs = list(nrow = control$bf.maxit)), 
            gz = replicate(p, matrix(nrow = n, ncol = control$bf.maxit), simplify = F)
        )
        names(trace.list$gz) <- xnames
    }
    while (rel.delta > control$bf.tol && cbf < control$bf.maxit){
        r <- y - beta0 - (gs %*% betas)
        deltaf <- 0
        for (j in 1:p){ #! /!\ Oscillation of betas and gz
            rj <- r + betas[j]*gs[,j]
            jterm <- single.term(x = x[[j]], y = rj, w = w, alpha = alphas[[j]], beta1 = betas[j], smoother = smoother[j], smoother.args = smoother.args[[j]], control = control)
            betas[j] <- jterm$beta 
            alphas[[j]] <- jterm$alpha
            deltaf <- deltaf + weighted.mean((jterm$gz - gs[,j])^2, w)  
            gs[,j] <- jterm$gz - mean(jterm$gz)
            beta0 <- beta0 + betas[j]*mean(jterm$gz)
        }
        rel.delta <- sqrt(deltaf/sum(apply(gs^2, 2, weighted.mean, w = w))) 
        cbf <- cbf + 1
        if (trace){ # Tracing the algorithm evolution
           trace.list$criteria[cbf,] <- c(deltaf, rel.delta, sum(r^2))
           trace.list$beta[cbf,] <- betas
           for (i in 1:p){
               trace.list$alpha[[i]][cbf,] <- alphas[[i]]
               trace.list$gz[[i]][,cbf] <- gs[,i]
           } 
        }
    }
    if (cbf == control$bf.maxit) warning(sprintf("Convergence not attained after %i iterations", control$bf.maxit))
    zs <- mapply("%*%", x, alphas)
    output <- list(alpha = alphas, gz = gs, z = zs, beta = c(beta0, betas))
    if (trace) output$trace <- trace.list
    return(output)
}

#! Generalized additive index models
aim <- function(x, y, w, family = gaussian, smoother = "spline", smoother.args = list(NULL), control = NULL, init.type = "runif", trace = T) 
# x: a list of matrices giving the variables for each index. The matrices need have the same number of rows (individuals) but can have diffeent numbers of columns (variables). If a matrix is provided, a single-index model is fitted. 
# y: a numeric vector containing the response observations
# w: weigths
# smoother: a character vector giving the smoother used for each index. Can be different for each index. Recycled if necessary.
# smoother.args: a (nested) list containing optional arguments for each smoother. 
# control: control parameters for the search algorithm
# trace: if true, keep criteria and parameters values of each iteration