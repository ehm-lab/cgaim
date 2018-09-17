#! The method of Yu and Ruppert (2002)
#! Difference: use of a LM instead of NLS for initial beta values
si.yr <- function(y, x, smoother = "ns", smoother.pars = list(df = 10), lambda = 1) 
{
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    # initialize alphas
    alpha <- alpha.init(type = "ols", p = p, y = y, x = x, normalize = TRUE, first.pos = TRUE)
    names(alpha) <- sprintf("alpha%i", 1:p)
    calpha <- alpha[-1] / alpha[1]
    z <- x %*% alpha
    smoother.pars$x <- z
    bases <- do.call(smoother, smoother.pars)
    nb <- ncol(bases)
    # augment data for penalized splines
    d <- diag(c(rep(0,1+p), rep(1,nb-p-1)))
    x.aug <- rbind(bases, sqrt(lambda)*d^(1/2))
    y.aug <- c(y, rep(0, nb))
    aug.data <- data.frame(y = y.aug, x = x.aug)
    # fit penalized splines for initial basis coefficients
    bcoef.init <- coef(lm(y ~ 0 + ., data = aug.data))
    # Define objective function
    obj.fun <- function(y, x, lambda, d, calpha, beta, smoother.pars){
        z <- x %*% param.map(calpha)
        smoother.pars$x <- z
        basis <- do.call(smoother, smoother.pars)
        nb <- ncol(basis)        
        pen <- sqrt(lambda)*d
        aug.dat <- rbind(basis, pen)
        aug.y <- c(y, rep(0, nb))
        aug.y - aug.dat%*%beta
    }
    # Fit NLS
    nls.fit <- nls(~ obj.fun(y, x, lambda, d, alpha, beta, smoother.pars), data = list(y = y, x = x, lambda = lambda, d = d, smoother.pars = smoother.pars), start = list(alpha = calpha, beta = bcoef.init), control = nls.control(warnOnly = T))  
    final.alpha <- param.map(coef(nls.fit)[1:(p-1)])
    final.z <- x %*% final.alpha
    final.beta <- coef(nls.fit)[p:length(coef(nls.fit))]
    smoother.pars$x <- final.z
    Bmat <- do.call(smoother, smoother.pars)
    final.gz <- Bmat %*% final.beta
    # GCV criterion (ou K-fold?)
#    cminv <- solve(t(Bmat)%*%Bmat + n*lambda*d)
#    As <- Bmat %*% cminv %*% t(Bmat)
#    trHat <- sum(diag(cminv%*%t(Bmat)%*%Bmat)) + sum(diag(solve(t(Bmat)%*%(diag(1,n)-As)%*%Bmat) %*% t(Bmat)%*%(diag(1,n) - As)^2%*%Bmat))
    output <- list(alpha = final.alpha, z = final.z, gz = final.gz, basis_coefs = final.beta)
    return(output)
}

param.map <- function(gamma){
    c(1, gamma) / sqrt(1 + sum(gamma^2))
}

p.spline.expansion <- function(expr, degree = 3, knots = NULL, nknots = 9, intercept = TRUE, coef.name = "beta")
# The one of Yu and Ruppert (2002)
# Same arguments as ns
{
    if (is.null(knots)){
       knots <- seq(0,1, length.out = nknots + 2)[-c(1,nknots+2)]
    }
    polys <- sprintf("(%s)^%i", expr, 1:degree)
    splines <- sprintf("((%1$s) - quantile(%1$s, %2$f))^%3$i", expr, knots, degree)
    basis <- c(polys, splines)
    if (intercept) basis <- c("1", basis)
    nb <- length(basis)
    coefs <- sprintf("%s%i", coef.name, 1:(nb+intercept)-intercept)
    cbprod <- paste(coefs, basis, sep = "*")
    final.form <- paste(cbprod, collapse = "+")
    return(list(formula = as.formula(final.form), coef.names = coefs))
}

# An additive version of Yu and Ruppert (2002) SI model
aim.yr <- function(x, y, w, smoother = "ns", smoother.args = list(NULL), smoother.penalty = list(NULL), lambdas = 1, control, init.type = "runif", trace = T)
# x: a list of matrices giving the variables for each index. The matrices need have the same number of rows (individuals) but can have diffeent numbers of columns (variables). If a matrix is provided, a single-index model is fitted. 
# y: a numeric vector containing the response observations
# w: weigths (TODO!)
# smoother: a character vector giving the smoother used for each index. Can be different for each index. Recycled if necessary.
# smoother.args: a (nested) list containing optional arguments for each smoother.
# smoother.penalty: a list containing penalty matrices for each index. 
# lambdas: a vector of penaly coefficients for each index
# control: control parameters for nls algorithm. See ?nls.control
# init.type: intialization of alphas TODO
# trace: if true, keep criteria and parameters values of each iteration TODO
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
    if (length(smoother.penalty) != p){
       warning("Smoother.penalty does not have the same length as x and is thus recycled.")
       smoother.penalty <- rep_len(smoother.penalty, p)
    }
    lambdas <- rep_len(lambdas, p)
    if(missing(w)) w <- rep(1/n,n)    
    # initialize alphas
    X <- Reduce(cbind, x)
    alpha <- alpha.init(type = "ols", p = p, y = y, x = X, normalize = TRUE, first.pos = F)
    names(alpha) <- NULL
    alphas <- split(alpha, rep(1:p, pvec))
    names(alphas) <- xnames
    calpha <- sapply(alphas, function(x) x[-1] / x[1])
    # Compute first z
    zs <- mapply("%*%", x, alphas, SIMPLIFY = F)
    # initialize betas
    bases <- vector("list", p)
    nb <- vector("integer", p)
    for (j in 1:p){
        smoother.args[[j]]$x <- zs[[j]]
        bases[[j]] <- do.call(smoother[j], smoother.args[[j]])
        nb[j] <- ncol(bases[[j]])
    }
    allbases <- Reduce(cbind, bases)
    for (j in 1:p) if(is.null(smoother.penalty[[j]])) smoother.penalty[[j]] <- diag(c(rep(0,1+pvec[j]), rep(1,nb[j]-pvec[j]-1)))
    Smat <- mapply("*", smoother.penalty, lambdas, SIMPLIFY = F)
    allS <- as.matrix(bdiag(Smat))
    aug.data <- data.frame(y = c(y, rep(0, sum(nb))), x = rbind(allbases, mat.sqrt(allS)))
    bcoef.init <- coef(lm(y ~ 0 + ., data = aug.data))
    # Define objective function
    obj.fun <- function(y, x, allS, calpha, beta, smoother.args){
        calpha <- split(calpha, rep(1:p, pvec-1))
        alphas <- lapply(calpha, param.map)
        zs <- mapply("%*%", x, alphas, SIMPLIFY = F)
        bases <- vector("list", p)
        for (j in 1:p){
            smoother.args[[j]]$x <- zs[[j]]
            bases[[j]] <- do.call(smoother[j], smoother.args[[j]])
        }
        allbases <- Reduce(cbind, bases)     
        aug.dat <- rbind(allbases, allS)
        aug.y <- c(y, rep(0, sum(nb)))
        aug.y - aug.dat%*%beta
    }
    # Fit NLS
    nls.fit <- nls(~ obj.fun(y, x, allS, calpha, beta, smoother.args), data = list(y = y, x = x, allS = allS, smoother.args = smoother.args), start = list(calpha = unlist(calpha), beta = bcoef.init), control = nls.control(warnOnly = T))
    fit.coef <- coef(nls.fit)
    final.calpha <- fit.coef[grep('calpha', names(fit.coef))]
    names(final.calpha) <- NULL
    final.calpha <- split(final.calpha, rep(1:p, pvec-1))
    final.alpha <- lapply(final.calpha, param.map)
    names(final.alpha) <- xnames
    final.z <- mapply("%*%", x, final.alpha)
    final.bases <- vector("list", p)
    for (j in 1:p){
        smoother.args[[j]]$x <- final.z[,j]
        final.bases[[j]] <- do.call(smoother[j], smoother.args[[j]])
    }
    final.beta <- split(fit.coef[grep('beta', names(fit.coef))], rep(1:p,nb))
    final.gz <- mapply("%*%", final.bases, final.beta)
    final.gz <- scale(final.gz, scale = F)
    beta0 <- sum(attr(final.gz, "scaled:center"))
    betas <- apply(final.gz, 2, function(x) sqrt(sum(x^2)))
    for (j in 1:p) final.gz[,j] <- final.gz[,j] / betas[j]
    output <- list(alpha = final.alpha, gz = final.gz, z = final.z, beta = c(beta0,betas))
    return(output)
}

#! Matrix square root
mat.sqrt<-function(S) # A simple matrix square root
{ 
    d<-eigen(S,symmetric=TRUE)
    rS<-d$vectors%*%diag(d$values^0.5)%*%t(d$vectors)
    return(rS)
}

#! A naïve AIM that attempts to integrate the remarks of Wang et al. (2015) that the index coefficients can be retrieved by LS
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

#! Trying to code myself the smoothing a Pya and Wood (2015). Fail
smooth.conspline <- function(x, y = NULL, nknots = 10, lambda = 1, ord = 4, knots, ...){
    n <- length(x)
    if (missing(knots)){ # conpute knots of the spline
       knots <- rep(0, ord + nknots + 2)
       knots[(ord + 2):(nknots + 1)] <- seq(min(x), max(x), length = nknots - ord)
       for (i in 1:(ord + 1)) knots[i] <- knots[ord + 2] - (ord + 2 - i) * (knots[ord + 3] - knots[ord + 2])
       for (i in (nknots + 2):(nknots + ord + 2)) knots[i] <- knots[nknots + 1] + (i - nknots - 1) * (knots[ord + 3] - knots[ord + 2])
    }
    # Prepare spline basis       
    B <- splineDesign(knots, x, ord = ord)
    Sigma <- matrix(1, nknots - 1, nknots - 1)
    Sigma[upper.tri(Sigma)] <- 0
    X <- B[,2:nknots] %*% Sigma
    nb <- ncol(X)
    P <- diff(diag(nknots - 1), difference = 1)
    P2 <- crossprod(P)
    M <- list(y = y, w = rep(1, n), X = X, C = matrix(0, 0, 0), S = list(P2), off = 0, sp = lambda, p = rep(0.1, nb), Ain = diag(nb), bin = c(-1e+12, rep(1e-12, nb - 1)))
    beta.tilde <- pcls(M)
    beta <- beta.tilde
    beta[-1] <- log(beta[-1])
}