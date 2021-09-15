
cgam_lookup <- c(inc = "s.incr", dec = "s.decr", cvx = "s.conv", ccv = "s.conc",
  inccvx = "s.incr.conv", deccvx = "s.decr.conv", incccv = "s.incr.conc",
  decccv = "s.decr.conc"
)


cgam.setup <- function(mf, control)
{
  # Prepare the formula
  mt <- attr(mf, "terms")
  ptot <- length(attr(mt, "term.labels"))
  gind <- attr(mt, "specials")$g
  sind <- attr(mt, "specials")$s
  # Default parameters for control
  defcontrol <- list()
  # Match control
  mc <- names(defcontrol) %in% names(control)
  control <- c(control, defcontrol[!mc])
  # Create call to "s" for each term
  term_s <- sapply(mf[c(gind, sind)], function(x){
    sfun <- cgam_lookup[attr(x, "fcons")]
    if (length(sfun) == 0) sfun <- "s"
    cpars <- list(str2lang(sfun), str2lang(attr(x, "label")))
    cpars <- c(cpars, attr(x, "s_opts"))
    cl <- as.call(cpars)
    deparse(cl)
  })
  # Add linear terms
  term_lin <- attr(mt, "term.labels")[-(c(gind, sind) - 1)]
  # Put together
  rhsform <- paste(c(term_s, term_lin), collapse = " + ")
  control$formula <- stats::reformulate(c(term_s, term_lin), all.vars(mt)[1])
  return(control)
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
  return(list(intercept = beta0, gz = gx, dgz = dgx, edf = sum(gfit$edf0)))
}