#!#############################################################################
#!
#!                        Main GAIM functions                                  
#!
#!#############################################################################

source("C:/Users/masselpl/Documents/Recherche/2017-2019 - Post-doc/Programmes/R/2 - Indices PPR/gaim/Secondary functions.R")

#! fit single term ppr
single.term <- function(x, y, w, alpha = alpha.init(ncol(x)), beta1 = 1, smoother = "spline", smoother.args = NULL, control = NULL, alpha.control = NULL) 
# y: vector of response values
# x: matrix of predictor values
# w: weights
# beta1: initial value of beta1 (possibly inherited from previous passes by the backfitting)
# smoother: the smoother used
# control: control parameters for the search algorithm
{
    def.control <- list(tol = 5e-3, max.iter = 20, min.step.len = 0.1, halving = T)
    control <- c(control,def.control[!names(def.control) %in% names(control)])
    def.alpha.control <- list(monotone = NULL, monotone.type = "QP", 
      normalize = "L2")
    alpha.control <- c(alpha.control,def.alpha.control[!names(def.alpha.control) %in% names(alpha.control)])
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
        delta <- gn_update(r, xf, w, monotone = alpha.control$monotone, 
          monotone.type = alpha.control$monotone.type)
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
    alpha <- normalize(alpha, alpha.control$normalize) 
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
aim.backfit <- function(x, y, w, smoother = "spline", smoother.args = list(NULL), control = NULL, init.type = "runif", trace = T) 
# x: a list of matrices giving the variables for each index. The matrices need have the same number of rows (individuals) but can have diffeent numbers of columns (variables). If a matrix is provided, a single-index model is fitted. 
# y: a numeric vector containing the response observations
# w: weigths
# smoother: a character vector giving the smoother used for each index. Can be different for each index. Recycled if necessary.
# smoother.args: a (nested) list containing optional arguments for each smoother. 
# control: control parameters for the search algorithm
# trace: if true, keep criteria and parameters values of each iteration
{
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
        f0norm <- sum(apply(gs^2, 2, weighted.mean, w = w))
        for (j in 1:p){ #! /!\ Oscillation of betas and gz
            rj <- r + betas[j]*gs[,j]
            jterm <- single.term(x = x[[j]], y = rj, w = w, alpha = alphas[[j]], beta1 = betas[j], smoother = smoother[j], smoother.args = smoother.args[[j]], control = control)
            betas[j] <- jterm$beta 
            alphas[[j]] <- jterm$alpha
            deltaf <- deltaf + weighted.mean((jterm$gz - gs[,j])^2, w)  
            gs[,j] <- jterm$gz - mean(jterm$gz)
            beta0 <- beta0 + betas[j]*mean(jterm$gz)
        }
        rel.delta <- sqrt(deltaf/f0norm) 
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
gaim <- function(x, y, w, family = gaussian, smoother = "spline", smoother.args = list(NULL), control = NULL, init.type = "runif", trace = T, backfit = T) 
# x: a list of matrices giving the variables for each index. The matrices need have the same number of rows (individuals) but can have diffeent numbers of columns (variables). If a matrix is provided, a single-index model is fitted. 
# y: a numeric vector containing the response observations
# w: weigths
# family: the distribution family of the response
# smoother: a character vector giving the smoother used for each index. Can be different for each index. Recycled if necessary.
# smoother.args: a (nested) list containing optional arguments for each smoother. 
# control: control parameters for the search algorithm
# trace: if true, keep criteria and parameters values of each iteration
# Backfit: probably temporary. Does a complete backfitting after each update or only one pass?
{
    if (!is.list(x)) x <- list(x)
    y <- as.vector(y)
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
    # initalization
    famObj <- do.call(family, list())
    beta0 <- famObj$linkfun(mean(y))
    betas <- rep(1,p)
    alphas <- mapply(alpha.init, p = pvec, x = x, MoreArgs = list(y = y, type = init.type))
    gs <- matrix(0, n, p)
    # Update
    rel.delta <- control$bf.tol + 1
    cbf <- 0
    if (trace){ # Tracing the algorithm evolution
        trace.list <- list(
            criteria = matrix(NA, nrow = control$bf.maxit, ncol = 2, dimnames = list(NULL, c("delta", "rel.delta"))), 
            beta = matrix(NA, nrow = control$bf.maxit, ncol = p, dimnames = list(NULL, xnames)), 
            alpha = mapply(matrix, ncol = pvec, MoreArgs = list(nrow = control$bf.maxit)), 
            gz = replicate(p, matrix(nrow = n, ncol = control$bf.maxit), simplify = F)
        )
        names(trace.list$gz) <- xnames
    }
    while (rel.delta > control$bf.tol && cbf < control$bf.maxit){
        # Compute adjusted variable and weights
        eta <- beta0 + (gs %*% betas)
        mu <- famObj$linkinv(eta)
        adj.y <- eta + (y - mu)/famObj$mu.eta(eta)
        adj.weights <- 1/(famObj$mu.eta(eta)*famObj$variance(mu))  #! How to integrate prior weights?
        deltaf <- 0
        f0norm <- sum(apply(gs^2, 2, weighted.mean, w = w))
        if (backfit){
            aim.res <- aim(x, adj.y, adj.weights, smoother = smoother, smoother.args = smoother.args, control =  control, init.type = init.type, trace = F)
            for (j in 1:p) deltaf <- deltaf + weighted.mean((aim.res$gz[,j] - gs[,j])^2, w)
            gs <- aim.res$gz
            beta0 <- aim.res$beta[1]
            betas <- aim.res$beta[-1]
            alphas <- aim.res$alpha
        } else {
            r <- adj.y - eta
            for (j in 1:p){ #! /!\ Oscillation of betas and gz
                rj <- r + betas[j]*gs[,j]
                jterm <- single.term(x = x[[j]], y = rj, w = w, alpha = alphas[[j]], beta1 = betas[j], smoother = smoother[j], smoother.args = smoother.args[[j]], control = control)
                betas[j] <- jterm$beta 
                alphas[[j]] <- jterm$alpha
                deltaf <- deltaf + weighted.mean((jterm$gz - gs[,j])^2, w)  
                gs[,j] <- jterm$gz - mean(jterm$gz)
                beta0 <- beta0 + betas[j]*mean(jterm$gz)
            }
        }
        rel.delta <- sqrt(deltaf/f0norm) 
        cbf <- cbf + 1
        if (trace){ # Tracing the algorithm evolution
           trace.list$criteria[cbf,] <- c(deltaf, rel.delta)
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

#! An alternative function

aim.GaussNewton <- function(x, y, w, bs = "tp", control = NULL, init.type = "runif")   #! ADD POSSIBLE ARGUMENTS FOR GAM AND SCAM
#! ADD BETA
{
    require(mgcv)
    require(scam)
    require(tsgam) #! /!\ only on github
    y <- as.vector(y)
    if (!is.list(x)) x <- list(x)
    if (is.null(names(x))) names(x) <- sprintf("V%i", 1:p)
    xnames <- names(x)
    Xall <- Reduce(cbind, x)
    n <- length(y)
    p <- length(x)
    pvec <- sapply(x, ncol)
    pind <- rep(1:p, pvec)
    if (any(sapply(x, nrow) != n)) stop("Matrices in x must have number of rows equal to the length of y")
    if(missing(w)) w <- rep(1/n,n)
    bs <- rep_len(bs, p) 
    def.control <- list(tol = 5e-3, max.iter = 20, min.step.len = 0.1, halving = T)
    control <- c(control,def.control[!names(def.control) %in% names(control)])
    # Initalization
    beta0 <- mean(y)
    y <- y - beta0
    alpha <- mapply(alpha.init, p = pvec, x = x, MoreArgs = list(y = y, type = init.type, first.pos = TRUE), SIMPLIFY =FALSE)  #! MAY CHANGE RESULT. STUDY THE IMPACT.
    z <- mapply("%*%", x, alpha) 
    form.rhs <- sprintf("s(%s, bs = '%s')", colnames(z), bs) #! ADD POSSIBILITY TO USE SCAM AND ARGUMENTS FOR s()
    form.gam <- sprintf("y ~ %s", paste(form.rhs, collapse = " + "))
    #! ADD LINEAR PREDICTORS    
    if (any(bs %in% c("mpi", "mpd", "cx", "cv", "micx", "micv", "mdcx", "mdcv"))){
       gfit <- scam(as.formula(form.gam), data = data.frame(y = y, z))
       dmod <- lapply(1:p, derivative.scam, object = gfit)
       dgz <- sapply(dmod, "[[", "d")
    } else {
       gfit <- gam(as.formula(form.gam), data = data.frame(y = y, z))  
       dmod <- fderiv(gfit, newdata = data.frame(z))  #! USE FUNCTION derivative.scam FOR SCAM
       dgz <- sapply(dmod$derivatives, "[[", "deriv")
    }    
    gz <- predict(gfit, type = "terms")
    yhat <- predict(gfit, type = "response")
    l2 <- L2(y, yhat, w)
    eps <- (var(y) - l2)/var(y)
    # Gauss-Newton search
    c1 <- 1
    while(eps > control$tol && c1 <= control$max.iter){
        r <- y - yhat
        dmat <- dgz[,pind]
        Vmat <- Xall*dmat #! POTENTIALLY RIDGE REGRESSION OR LASSO? SEE ROOSEN AND HASTIE (1994) AND SEARCH LITTERATURE
        delta <- coef(lm(r ~ 0 + Vmat))
        cuts <- 1
        l2.old <- l2
        repeat {   #Loop for halving steps in case of bad step  #
            delta <- delta * cuts
            alpha.new <- mapply("+", alpha, split(delta, pind), SIMPLIFY = FALSE)
            z <- mapply("%*%", x, alpha.new)
            if (any(bs %in% c("mpi", "mpd", "cx", "cv", "micx", "micv", "mdcx", "mdcv"))){
               gfit <- scam(as.formula(form.gam), data = data.frame(y = y, z))
               dmod <- lapply(1:p, derivative.scam, object = gfit)
               dgz <- sapply(dmod, "[[", "d")
            } else {
               gfit <- gam(as.formula(form.gam), data = data.frame(y = y, z))  
               dmod <- fderiv(gfit, newdata = data.frame(z))  #! USE FUNCTION derivative.scam FOR SCAM
               dgz <- sapply(dmod$derivatives, "[[", "deriv")
            }    
            gz <- predict(gfit, type = "terms")
            yhat <- predict(gfit, type = "response")
            l2 <- L2(as.vector(y), yhat, w)            
            if (l2 < l2.old || cuts < control$min.step.len || !control$halving){
               break
            }
            cuts <- cuts / 2
        }
        alpha <- alpha.new
        eps <- (l2.old - l2)/l2.old
        c1 <- c1 + 1
    }
    alpha <- lapply(alpha, function(x) sign(x[1])*x/sqrt(sum(x^2))) #! FIRS ELEMENT POSITIVE? EACH ITERATION?
    z <- mapply("%*%", x, alpha)
    if (any(bs %in% c("mpi", "mpd", "cx", "cv", "micx", "micv", "mdcx", "mdcv"))){
       gfit <- scam(as.formula(form.gam), data = data.frame(y = y, z))
    } else {
       gfit <- gam(as.formula(form.gam), data = data.frame(y = y, z))  
    }    
    gz <- predict(gfit, type = "terms")
    beta1 <- apply(gz, 2, function(x) sqrt(sum(x^2)))
    gz <- gz / matrix(beta1, ncol = p, nrow = n, byrow = T)
    return(list(alpha = alpha, z = z, gz = gz, beta = c(beta0,beta1), fitted.values = beta0 + predict(gfit, type = "response")))
}




#! At each iteration of the Newton-Raphson, fit a whole AIM model
gaim2 <- function(x, y, w, family = gaussian, control = NULL, init.type = "runif", trace = T) 
# x: a list of matrices giving the variables for each index. The matrices need have the same number of rows (individuals) but can have diffeent numbers of columns (variables). If a matrix is provided, a single-index model is fitted. 
# y: a numeric vector containing the response observations
# w: weigths
# family: the distribution family of the response
# smoother: a character vector giving the smoother used for each index. Can be different for each index. Recycled if necessary.
# smoother.args: a (nested) list containing optional arguments for each smoother. 
# control: control parameters for the search algorithm
# trace: if true, keep criteria and parameters values of each iteration
# Backfit: probably temporary. Does a complete backfitting after each update or only one pass?
{
    require(mgcv)
    require(scam)
    require(tsgam) #! /!\ only on github
    y <- as.vector(y)
    if (!is.list(x)) x <- list(x)
    if (is.null(names(x))) names(x) <- sprintf("V%i", 1:p)
    xnames <- names(x)
    Xall <- Reduce(cbind, x)
    n <- length(y)
    p <- length(x)
    pvec <- sapply(x, ncol)
    pind <- rep(1:p, pvec)
    if (any(sapply(x, nrow) != n)) stop("Matrices in x must have number of rows equal to the length of y")
    if(missing(w)) w <- rep(1/n,n)
    def.control <- list(tol = 5e-3, max.iter = 50, min.step.len = 0.1, halving = T)
    control <- c(control,def.control[!names(def.control) %in% names(control)])
    # initalization
    famObj <- do.call(family, list())
    beta0 <- famObj$linkfun(mean(y))
    betas <- rep(1,p)
    alphas <- mapply(alpha.init, p = pvec, x = x, MoreArgs = list(y = y, type = init.type))
    gs <- matrix(0, n, p)
    # Update
    rel.delta <- control$tol + 1
    cbf <- 0
    while (rel.delta > control$tol && cbf < control$max.iter){
        # Compute adjusted variable and weights
        eta <- beta0 + (gs %*% betas)
        mu <- famObj$linkinv(eta)
        adj.y <- eta + (y - mu)/famObj$mu.eta(eta)
        adj.weights <- 1/(famObj$mu.eta(eta)*famObj$variance(mu))  #! How to integrate prior weights?
        deltaf <- 0
        f0norm <- sum(apply(gs^2, 2, weighted.mean, w = w))
        aim.res <- aim2(x, adj.y, adj.weights, control =  control, init.type = init.type)
        for (j in 1:p) deltaf <- deltaf + weighted.mean((aim.res$gz[,j] - gs[,j])^2, w)
        gs <- aim.res$gz
        beta0 <- aim.res$beta[1]
        betas <- aim.res$beta[-1]
        alphas <- aim.res$alpha
        rel.delta <- sqrt(deltaf/f0norm) 
        cbf <- cbf + 1
    }
    if (cbf == control$max.iter) warning(sprintf("Convergence not attained after %i iterations", control$bf.maxit))
    zs <- mapply("%*%", x, alphas)
    output <- list(alpha = alphas, gz = gs, z = zs, beta = c(beta0, betas))
    return(output)
}


#! At each iteration of the Newton-Raphson, make a single step of the Gauss-Newton algorithm
#! ÃƒÆ’Ã¢â€šÂ¬ repenser
gaim3 <- function(x, y, w, control = NULL, init.type = "runif")
{
    require(mgcv)
    require(scam)
    require(tsgam) #! /!\ only on github
    y <- as.vector(y)
    if (!is.list(x)) x <- list(x)
    if (is.null(names(x))) names(x) <- sprintf("V%i", 1:p)
    xnames <- names(x)
    Xall <- Reduce(cbind, x)
    n <- length(y)
    p <- length(x)
    pvec <- sapply(x, ncol)
    pind <- rep(1:p, pvec)
    if (any(sapply(x, nrow) != n)) stop("Matrices in x must have number of rows equal to the length of y")
    if(missing(w)) w <- rep(1/n,n)
    def.control <- list(tol = 5e-3, max.iter = 20, min.step.len = 0.1, halving = T)
    control <- c(control,def.control[!names(def.control) %in% names(control)])
    # Initalization
    famObj <- do.call(family, list())
    beta0 <- famObj$linkfun(mean(y))
    alpha <- mapply(alpha.init, p = pvec, x = x, MoreArgs = list(y = y, type = init.type, first.pos = TRUE))
    z <- mapply("%*%", x, alpha) 
    form.rhs <- sprintf("s(%s)", colnames(z)) #! ADD POSSIBILITY TO USE SCAM AND ARGUMENTS FOR s()
    form.gam <- sprintf("y ~ %s", paste(form.rhs, collapse = " + "))
    gfit <- gam(as.formula(form.gam), data = data.frame(y = y, z), family = family)  #! ADD LINEAR PREDICTORS
    gz <- predict(gfit, type = "terms")
    dmod <- fderiv(gfit, newdata = data.frame(z))  #! USE FUNCTION derivative.scam FOR SCAM
    dgz <- sapply(dmod$derivatives, "[[", "deriv")
    yhat <- predict(gfit, type = "response")
    l2 <- L2(y, yhat, w)
    eps <- (var(y) - l2)/var(y)
    # Gauss-Newton search
    c1 <- 1
    while(eps > control$tol && c1 <= control$max.iter){
        eta <- beta0 + (gs %*% betas)
        mu <- famObj$linkinv(eta)
        adj.y <- eta + (y - mu)/famObj$mu.eta(eta)
        adj.weights <- 1/(famObj$mu.eta(eta)*famObj$variance(mu))
        r <- adj.y - as.vector(predict(gfit, type = "link"))
        dmat <- dgz[,pind]
        Vmat <- Xall*dmat
        delta <- coef(lm(r ~ 0 + Vmat))
        cuts <- 1
        l2.old <- l2
        repeat {   #Loop for halving steps in case of bad step  #
            delta <- delta * cuts
            alpha.new <- mapply("+", alpha, split(delta, pind))
            z <- mapply("%*%", x, alpha.new)
            gfit <- gam(as.formula(form.gam), data = data.frame(y = y, z))#! ADD POSSIBILITY TO USE SCAM AND ARGUMENTS FOR s()
            gz <- predict(gfit, type = "terms")
            dmod <- fderiv(gfit, newdata = data.frame(z))  #! USE FUNCTION derivative.scam FOR SCAM
            dgz <- sapply(dmod$derivatives, "[[", "deriv")
            yhat <- predict(gfit, type = "response")
            l2 <- L2(as.vector(y), yhat, w)            
            if (l2 < l2.old || cuts < control$min.step.len || !control$halving){
               break
            }
            cuts <- cuts / 2
        }
        alpha <- alpha.new
        eps <- (l2.old - l2)/l2.old
        c1 <- c1 + 1
    }
    alpha <- sapply(alpha, function(x) sign(x[1])*x/sqrt(sum(x^2))) #! FIRS ELEMENT POSITIVE?
    z <- mapply("%*%", x, alpha)
    gfit <- gam(as.formula(form.gam), data = data.frame(y = y, z))#! ADD POSSIBILITY TO USE SCAM AND ARGUMENTS FOR s()
    gz <- predict(gfit, type = "terms")
    beta1 <- apply(gz, 2, function(x) sqrt(sum(x^2)))
    gz <- gz / matrix(beta1, ncol = p, nrow = n, byrow = T)
    return(list(alpha = alpha, z = z, gz = gz, beta = c(beta0,beta1)))
}