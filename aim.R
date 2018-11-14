#! A naïve AIM that attempts to integrate the remarks of Wang et al. (2015) that the index coefficients can be retrieved by LS
#! Steps:
#!    1. Estimate coefficients alpha through OLS of all X on Y
#!    2. Estimate ridge functions by GAM on the constructed indices
aim.naive <- function(x, y, w, gam.pars = list(), control, trace = T)
{
    if (!is.list(x)) x <- list(x)
    xnames <- names(x)
    n <- length(y)
    p <- length(x)
    pvec <- sapply(x, ncol)
    if (any(sapply(x, nrow) != n)) stop("Matrices in x must have number of rows equal to the length of y")
    if(missing(w)) w <- rep(1/n,n)    
    # estimate alphas
    X <- Reduce(cbind, x)
    alpha <- alpha.init(type = "ols", p = p, y = y, x = X, normalize = F, first.pos = F)
    names(alpha) <- NULL
    alphas <- split(alpha, rep(1:p, pvec))
    names(alphas) <- xnames
    nor.alphas <- lapply(alphas, function(x) x/sqrt(sum(x^2)))
    # compute zs
    zs <- mapply("%*%", x, nor.alphas)
    # estimate gs
    gam.pars$data <- data.frame(y = y, zs)
    gam.pars$formula <- as.formula(sprintf("y ~ %s", paste(sprintf("s(%s)", names(gam.pars$data)[-1]), collapse = "+")))
    gz.fit <- do.call(gam, gam.pars)
    gz <- predict(gz.fit, type = "terms")
    final.gz <- scale(gz, scale = F)
    beta0 <- sum(attr(final.gz, "scaled:center"))
    betas <- apply(final.gz, 2, function(x) sqrt(sum(x^2))) 
    for (j in 1:p) final.gz[,j] <- final.gz[,j] / betas[j]
    output <- list(alpha = alphas, gz = final.gz, z = zs, beta = c(beta0,betas))
    return(output)   
}


#! AIM through classical backfitting
#! Loop through all indices X1, ..., Xp, X1, ..., Xp, ...
#!    Estimate ridge function and coefficients alpha by single-index model.
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
    alphas <- mapply(alpha.init, p = pvec, x = x, MoreArgs = list(y = y, type = init.type), SIMPLIFY = F)
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

#! AIM by Gauss-Newton algorithm
#! Repeat two steps until convergence:
#!    1) Update coefficient of all indices at once by Gauss-Newton
#!    2) Apply GAM (or SCAM) on indices to update Ridge Functions.
aim.GaussNewton <- function(x, y, w, bs = "tp", control = NULL, init.type = "runif")   #! ADD POSSIBLE ARGUMENTS FOR GAM AND SCAM
#! ADD BETA
{
    require(mgcv)
    require(scam)
    require(gratia) #! /!\ only on github
    #! Remove the requirements
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
       dmod <- gratia::fderiv(gfit, newdata = data.frame(z))  #! github package gratia. May want to change this line if the package is not available on CRAN when the paper is published
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
        #! Replace by compute.update
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