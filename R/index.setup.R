# Sets up alpha estimation
index.setup <- function(mf, alpha_control)
{
  # Initialize controls
  alpha_control <- do.call(alpha.control, alpha_control)
  # Extract index variables
  mt <- attr(mf, "terms")
  index_interp <- mf[attr(mt, "specials")$g]
  # Labels for smoothing step
  index_labels <- make.names(sapply(index_interp, attr, "label"), unique = T)
  # Extract response
  y <- stats::model.response(mf)
  attr(y, "varname") <- all.vars(mt)[1]
  # Extract weights
  w <- mf$`(weights)`
  if (is.null(w)) w <- rep(1, NROW(y))
  # Organizing variables in indices as a matrix with index positions   
  p <- length(index_interp)
  pvec <- sapply(index_interp, attr, "nterms")
  index <- rep(1:p, pvec)
  Xind <- do.call(cbind, index_interp)
  names(index) <- rep(index_labels, pvec)
  # Number of variables and indices
  ptot <- sum(sapply(index_interp, ncol))
  # Check Cmat
  if (!is.null(alpha_control$Cmat)){
    if (is.vector(alpha_control$Cmat)) alpha_control$Cmat <- 
        t(as.matrix(alpha_control$Cmat))
    if (ncol(alpha_control$Cmat) != ptot){
      stop("Number of columns in alpha_control$Cmat does not match
        the indices")
    }
  }
  # Add Cmat provided in index specific specifications
  alpha_control$Cmat <- rbind(alpha_control$Cmat, 
    as.matrix(Matrix::bdiag(lapply(index_interp, attr, "Cmat"))))
  # Add bvec
  alpha_control$bvec <- c(alpha_control$bvec, 
    unlist(lapply(index_interp, attr, "bvec"), use.names = F))
  # Check if there are any duplicates
  dups <- duplicated(alpha_control$Cmat)
  alpha_control$Cmat <- alpha_control$Cmat[!dups,,drop = F]
  alpha_control$bvec <- alpha_control$bvec[!dups]
  # Check number of constraints
  dims <- dim(alpha_control$Cmat)
  if (dims[1] > dims[2]) warning("More constraints than parameters. Check
    possible redundant constraints.")
  # Initialize alpha
  if (is.null(alpha_control$alpha.start)){
    init_pars <- alpha_control
    init_pars <- c(init_pars, list(y = y, x = Xind, w = w, index = index))
    alpha_control$alpha <- do.call(alpha_init, init_pars)
  } else {
    alpha_control$alpha <- unlist(alpha_control$alpha.start)
    if(length(alpha_control$alpha) != ptot){
      warning(paste0("alpha.start length is inconsistent with index matrix ", 
        "and is recycled"))
      alpha_control$alpha <- rep_len(alpha_control$alpha, ptot)
    }
  }
  alpha_control$alpha.start <- alpha_control$init.type <- NULL
  return(list(y = y, x = Xind, index = index, w = w,
    alpha_control = alpha_control)) 
}