
scam.setup <- function(mt, data, smooth.control){
  # Prepare the covariates and formula
  ptot <- length(attr(mt, "term.labels"))
  gind <- attr(mt, "specials")$g
  smooth_terms <- c(all.vars(mt)[1], vector("character", ptot))
  index_interp <- get("index_interp", envir = parent.frame())
  index_forms <- sapply(index_interp, attr, "opt_formula")
  smooth_terms[gind] <- sprintf("s(%s)", 
    mapply(paste, get("index_labels", envir = parent.frame()), 
      index_forms, MoreArgs = list(sep = ", ")))
  smooth_terms[-c(1, gind)] <- attr(mt[-(gind - 1)], "term.labels")
  gam_formula <- stats::reformulate(smooth_terms[-1],
     response = smooth_terms[1])
  smooth.control$formula <- gam_formula
  covariate_names <- setdiff(all.vars(mt)[-1], 
    unlist(lapply(index_interp, attr, "term")))
  smooth.control$Xcov <- data[,covariate_names]
  # Special default parameters (tol)
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
  gx <- gx[,match(colnames(scam_data)[1:p], colnames(gx))]
  segx <- segx[,match(colnames(scam_data)[1:p], colnames(segx))]
  # Estimate first derivative (different functions for scam or gam)
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

scar.setup <- function(mt, data, smooth.control){
  smooth.control <- smooth.control[names(smooth.control) %in% 
    methods::formalArgs(scar::scar)]
  gind <- attr(mt, "specials")$g
  ptot <- length(attr(mt, "term.labels"))
  covind <- (1:ptot)[-(gind - 1)]
  if (length(covind) > 0){
    mtcov <- mt[-(gind - 1)]
    mfcov <- stats::model.frame(
      stats::reformulate(all.vars(mtcov)[-1], intercept = F), 
      data = data)
    if (!all(sapply(mfcov,is.numeric))){
      stop("All variables used in 'scar' must be numeric")
    }
    smooth.control$Xcov <- do.call(get("na.action", envir = parent.frame()), 
      list(object = stats::model.matrix(mfcov, data = data)))
  }
  bss <- sapply(1:ptot, function(i){
    as.list(mt[i][[3]])$bs
  })
  notbs <- which(sapply(bss, is.null))
  sind <- attr(mt, "specials")$s
  if (length(intersect(notbs, c(gind - 1, sind - 1))) > 0){
    stop("All nonlinear terms must be shape-constrained when using 'scar'")
  } 
  bss[notbs] <- "l"
  smooth.control$shape <- c(unlist(bss[gind - 1]), unlist(bss[covind]))
  return(smooth.control)
}

smooth_scar <- function(x, y, shape, Xcov = NULL, ...)
{
  scar_data <- cbind(x, Xcov)
  gfit <- scar::scar(x = scar_data, y = y, shape = shape, ...)
  gx <- gfit$componentfit
  dgx <- mapply(function(x, gx) stats::splinefun(x, gx)(x, deriv = 1), 
    as.data.frame(scar_data), as.data.frame(gx))
  beta0 <- gfit$constant
  return(list(intercept = beta0, gz = gx, dgz = dgx))
}

