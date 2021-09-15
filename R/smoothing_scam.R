scam_lookup <- c(inc = "mpi", dec = "mpd", cvx = "cx", ccv = "cv",
  inccvx = "micx", deccvx = "mdcx", incccv = "micv",
  decccv = "mdcv"
)

scam.setup <- function(mf, control)
{
  # Prepare the formula
  mt <- attr(mf, "terms")
  ptot <- length(attr(mt, "term.labels"))
  gind <- attr(mt, "specials")$g
  sind <- attr(mt, "specials")$s
  # Default parameters for control
  fcons <- sapply(mf[c(gind, sind)], "attr", "fcons")
  defcontrol <- list()
  if (length(unlist(fcons)) > 0) defcontrol$sp <- rep(0, length(c(gind, sind)))
  # Match control
  mc <- names(defcontrol) %in% names(control)
  control <- c(control, defcontrol[!mc])
  # Create call to "s" for each term
  term_s <- sapply(mf[c(gind, sind)], function(x){
    cpars <- list(str2lang("s"), str2lang(attr(x, "label")))
    fcons <- attr(x, "fcons")
    if (!is.null(fcons)) cpars$bs <- unname(scam_lookup[fcons])
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


smooth_scam <- function(x, y, formula, Xcov, ...)
{ 
  # Construct data for scam fitting
  scam_data <- data.frame(as.data.frame(x), Xcov)
  p <- ncol(scam_data)
  n <- nrow(scam_data)
  scam_data[attr(y, "varname")] <- y
  iscons <- any(sapply(scam_lookup, grepl, as.character(formula)[3]))
  if (iscons){
    gfit <- scam::scam(formula, data = scam_data, ...)
    # print("scam"); flush.console()
  } else {
    gfit <- mgcv::gam(formula, data = scam_data, ...)
    derivs <- gratia::fderiv(gfit, newdata = scam_data)
    # print("gam"); flush.console()
  }
  # Extract estimated terms
  gx <- stats::predict(gfit, type = "terms")
  colnames(gx) <- gsub("s\\(|\\)", "", colnames(gx))
  gx <- gx[,match(colnames(scam_data)[1:p], colnames(gx)), drop = FALSE]
  # Estimate first derivative
  trms <- stats::terms(formula, specials = "s") 
  smterms <- attr(trms, "specials")$s - 1
  dgx <- matrix(0, n, p)
  for (j in 1:p){
    jind <- which(all.vars(formula)[-1] == names(scam_data)[j])
    if (jind %in% smterms){
      if (iscons){
        dgx[,j] <- scam::derivative.scam(gfit, jind)$d
      } else {
        dgx[,j] <- derivs$derivatives[[which(smterms == jind)]]$deriv
      }
    } else {
      if (is.numeric(scam_data[j])){
        dgx[,j] <- stats::coef(gfit)[names(scam_data)[j]]
      }
    }
  }
  beta0 <- stats::coef(gfit)[1]
  return(list(intercept = beta0, gz = gx, dgz = dgx, edf = sum(gfit$edf)))
}
