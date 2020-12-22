#-----------------------------------
# cgaim specific functions
#-----------------------------------

bspl <- function(x, knots = NULL, df = 10, ord = 4, intercept = FALSE)
{
  rng <- range(x)
  if (is.null(knots)){
    nk <- df - ord + 1
    dk <- diff(rng) / nk
    knots <- seq(rng[1] - dk * (ord - 1), rng[2] + dk * (ord - 1), by = dk)
  }
  bsmat <- splines::splineDesign(knots = knots, x = x, ord = ord, 
    outer.ok = TRUE)
  if (!intercept){
    bsmat <- bsmat[,-1]
  }
  return(bsmat)
}


#-----------------------------------
# scam package
#-----------------------------------

scam_lookup <- c(inc = "mpi", dec = "mpd", cvx = "cx", ccv = "cv",
  inccvx = "micx", deccvx = "mdcx", incccv = "micv",
  decccv = "mdcv"
)

scam.setup <- function(mt, data, smooth.control){
  # Check parameters
  scamargs <- methods::formalArgs(scam::scam)
  if (length(setdiff(names(smooth.control), scamargs)) > 0){
    warning("'smooth.control' contains arguments unused by the 'scam' function")
  }
  smooth.control <- smooth.control[names(smooth.control) %in% scamargs]
  # Prepare the covariates and formula
  ptot <- length(attr(mt, "term.labels"))
  gind <- attr(mt, "specials")$g
  form_terms <- c(all.vars(mt)[1], attr(mt, "term.labels"))
  index_interp <- stats::model.frame(mt[gind - 1], data = data)[-1]
  index_forms <- sapply(index_interp, attr, "label")
  index_fcons <- lapply(index_interp, attr, "fcons")
  whichfcons <- !sapply(index_fcons, is.null) 
  index_forms[whichfcons] <- sprintf("%s, bs = '%s'", index_forms[whichfcons], 
    scam_lookup[unlist(index_fcons[whichfcons])])
  index_pars <- sapply(index_interp, attr, "s_formula")
  whichpars <- nchar(index_pars) > 0
  index_forms[whichpars] <- sprintf("%s, %s", index_forms[whichpars], 
    index_pars[whichpars])
  form_terms[gind] <- sprintf("s(%s)", index_forms)
  # Covariate formula
  sind <- attr(mt, "specials")$s
  sp <- length(sind)
  if (sp > 0){
    sinterp <- stats::model.frame(mt[sind - 1], data = data)[-1]
    slabs <- sapply(sinterp, attr, "label")
    scons <- lapply(sinterp, attr, "fcons")
    snonull <- !sapply(scons, is.null)
    slabs[snonull] <- sprintf("%s, bs = '%s'", slabs[snonull], 
      scam_lookup[unlist(scons[snonull])])
    sopts <- sapply(sinterp, attr, "s_formula")
    slabs[nchar(sopts) > 0] <- paste(slabs[nchar(sopts) > 0],
      sopts[nchar(sopts) > 0], sep = ", ")
    sform <- sprintf("s(%s)", slabs)
    form_terms[sind] <- sform
  }
  # final object
  smooth.control$formula <- stats::reformulate(form_terms[-1], form_terms[1])
  covar <- all.vars(mt[-(gind - 1)])[-1]
  smooth.control$Xcov <- data[covar]
  return(smooth.control)
}

smooth_scam <- function(x, y, formula, Xcov, ...)
{ 
  # Construct data for scam fitting
  scam_data <- data.frame(as.data.frame(x), Xcov)
  p <- ncol(scam_data)
  n <- nrow(scam_data)
  scam_data[attr(y, "varname")] <- y
  gfit <- scam::scam(formula, data = scam_data, ...)
  # Extract estimated terms
  gxp <- stats::predict(gfit, type = "terms", se.fit = T)
  gx <- gxp$fit
  segx <- gxp$se.fit
  colnames(gx) <- colnames(segx) <- gsub("s\\(|\\)", "", colnames(gx))
  gx <- gx[,match(colnames(scam_data)[1:p], colnames(gx)), drop = FALSE]
  segx <- segx[,match(colnames(scam_data)[1:p], colnames(segx)), drop = FALSE]
  # Estimate first derivative 
  smterms <- attr(stats::terms(formula, specials = "s"),
    "specials")$s - 1
  dgx <- matrix(0, n, p)
  for (j in 1:p){
    jind <- which(all.vars(formula)[-1] == names(scam_data)[j])
    if (jind %in% smterms){
      dgx[,j] <- scam::derivative.scam(gfit, jind)$d
    } else {
      if (is.numeric(scam_data[j])){
        dgx[,j] <- stats::coef(gfit)[names(scam_data)[j]]
      }
    }
  }
  beta0 <- stats::coef(gfit)[1]
  return(list(intercept = beta0, gz = gx, dgz = dgx, se = segx))
}

#-----------------------------------
# scar package
#-----------------------------------

scar_lookup <- c(inc = "in", dec = "de", cvx = "cvx", ccv = "ccv",
  inccvx = "cvxin", deccvx = "cvxde", incccv = "ccvin",
  decccv = "ccvde"
)

scar.setup <- function(mt, data, smooth.control){
  smooth.control <- smooth.control[names(smooth.control) %in% 
    methods::formalArgs(scar::scar)]
  # Indices
  gind <- attr(mt, "specials")$g
  ptot <- length(attr(mt, "term.labels"))
  gsh <- rep("l", length(gind))
  index_interp <- stats::model.frame(mt[gind - 1], data = data)[-1]
  fcons <- lapply(index_interp, attr, "fcons")
  nonull <- !sapply(fcons, is.null)
  gsh[nonull] <- scar_lookup[unlist(fcons[nonull])]
  smooth.control$shape <- gsh
  lflag <- any(gsh == "l")
  # Covariates
  covind <- (1:ptot)[-(gind - 1)]
  if (length(covind) > 0){
    mtcov <- mt[-(gind - 1)]
    mfcov <- stats::model.frame(
      stats::reformulate(all.vars(mtcov)[-1], intercept = F), 
      data = data)
    if (!all(sapply(mfcov,is.numeric))){
      stop("All variables used in 'scar' must be numeric")
    }
    smooth.control$Xcov <- mfcov
    # Covariate shape
    ssh <- rep("l", ncol(mfcov))
    sind <- attr(mtcov, "specials")$s - 1
    sinterp <- stats::model.frame(mtcov[sind], data = data)[-1]
    scons <- lapply(sinterp, attr, "fcons")
    snonull <- !sapply(scons, is.null)
    ssh[sind[snonull]] <- scar_lookup[unlist(scons[snonull])]
    smooth.control$shape <- c(smooth.control$shape, ssh)
    lflag <- lflag || (length(intersect(sind, which(ssh == "l"))) > 0 )
  }
  if (lflag){
    warning("All nonlinear terms must be shape-constrained when using 'scar'. Unconstrained smooth terms have been set as linear.")
  } 
  return(smooth.control)
}

smooth_scar <- function(x, y, shape, Xcov = NULL, ...)
{
  scar_data <- data.matrix(cbind(x, Xcov))
  gfit <- scar::scar(x = scar_data, y = y, shape = shape, ...)
  gx <- gfit$componentfit
  dgx <- mapply(function(x, gx) stats::splinefun(x, gx)(x, deriv = 1), 
    as.data.frame(scar_data), as.data.frame(gx))
  beta0 <- gfit$constant
  return(list(intercept = beta0, gz = gx, dgz = dgx))
}

#-----------------------------------
# cgam package
#-----------------------------------

cgam_lookup <- c(inc = "s.incr", dec = "s.decr", cvx = "s.conv", ccv = "s.conc",
  inccvx = "s.incr.conv", deccvx = "s.decr.conv", incccv = "s.incr.conc",
  decccv = "s.decr.conc"
)

cgam.setup <- function(mt, data, smooth.control){
  cgamargs <- methods::formalArgs(cgam::cgam)
  if (length(setdiff(names(smooth.control), cgamargs)) > 0){
    warning("'smooth.control' contains arguments unused by the 'cgam' function")
  }
  smooth.control <- smooth.control[names(smooth.control) %in% cgamargs]
  form_terms <- c(all.vars(mt)[1], attr(mt, "term.labels"))
  # g terms
  gind <- attr(mt, "specials")$g
  gp <- length(gind)
  index_interp <- stats::model.frame(mt[gind - 1], data = data)[-1]
  gcons <- lapply(index_interp, attr, "fcons")
  gnonull <- !sapply(gcons, is.null)
  gfuns <- rep("s", gp)
  gfuns[gnonull] <- cgam_lookup[unlist(gcons[gnonull])]
  gopts <- sapply(index_interp, attr, "s_formula")
  gform <- sapply(index_interp, attr, "label")
  gform[nchar(gopts) > 0] <- paste(gform[nchar(gopts) > 0],
    gopts[nchar(gopts) > 0], sep = ", ")
  gform <- sprintf("%s(%s)", gfuns, gform)
  form_terms[gind] <- gform
  # s terms
  sind <- attr(mt, "specials")$s
  sp <- length(sind)
  if (sp > 0){
    sinterp <- stats::model.frame(mt[sind - 1], data = data)[-1]
    scons <- lapply(sinterp, attr, "fcons")
    snonull <- !sapply(scons, is.null)
    sfuns <- rep("s", sp)
    sfuns[snonull] <- cgam_lookup[unlist(scons[snonull])]
    sopts <- sapply(sinterp, attr, "s_formula")
    sform <- sapply(sinterp, attr, "label")
    sform[nchar(sopts) > 0] <- paste(sform[nchar(sopts) > 0],
      sopts[nchar(sopts) > 0], sep = ", ")
    sform <- sprintf("%s(%s)", sfuns, sform)
    form_terms[sind] <- sform
  }
  # final object
  smooth.control$formula <- stats::reformulate(form_terms[-1], form_terms[1])
  covar <- all.vars(mt[-(gind - 1)])[-1]
  smooth.control$Xcov <- data[covar]
  return(smooth.control)
}

smooth_cgam <- function(x, y, formula, Xcov, ...){
  cgam_data <- data.frame(as.data.frame(x), Xcov)  
  p <- ncol(cgam_data)
  n <- nrow(cgam_data)
  cgam_data[attr(y, "varname")] <- y
  gfit <- cgam::cgam(formula, data = cgam_data, ...)
  # Extract estimated terms
  preds <- t(gfit$etacomp)
  if (!is.null(gfit$zmat)){ 
    preds <- cbind(preds, mapply("*", as.data.frame(gfit$zmat), gfit$zcoef[-1]))
  }
  colnames(preds) <- c(gfit$xnms, gfit$znms)
  whichinds <- pmatch(colnames(x), colnames(preds))
  gx <- preds[,c(whichinds, setdiff(seq_len(ncol(preds)), whichinds)), drop = FALSE]
  # Estimate first derivative (different functions for scam or gam)
  xmat <- cbind(gfit$xmat_add, gfit$zmat)
  colnames(xmat) <- c(gfit$xnms, gfit$znms)
  xmat <- xmat[,c(whichinds, setdiff(seq_len(ncol(preds)), whichinds)), drop = FALSE]
  dgx <- mapply(function(x, gx) stats::splinefun(x, gx)(x, deriv = 1), 
    as.data.frame(gfit$xmat_add), as.data.frame(t(gfit$etacomp)))
  nz <- length(gfit$zcoef[-1])
  if (nz > 0) dgx <- cbind(dgx, matrix(gfit$zcoef[-1], nrow = n, ncol = nz))
  colnames(dgx) <- c(gfit$xnms, gfit$znms)
  dgx <- dgx[,c(whichinds, setdiff(seq_len(ncol(preds)), whichinds)), drop = FALSE]
  beta0 <- gfit$coef[1]
  return(list(intercept = beta0, gz = gx, dgz = dgx, se = NULL))
}