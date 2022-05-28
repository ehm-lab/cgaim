#' Plot CGAIM terms
#'
#' Plot method for the ridge and smooth terms of a \code{cgaim} object. If
#'    provided, may also plot confidence intervals.
#'
#' At the moment, only plots ridge and smooth functions. When several terms
#'    are plotted, waits for user response before displaying the next one.
#'
#' @param x A \code{cgaim} object.
#' @param select A numeric or character vector indicating which terms
#'    to plot.
#' @param ci An object returned by a call to \code{\link{confint.cgaim}}. If
#'    \code{NULL}, no confidence interval is drawn.
#' @param ci.plot Whether to plot the confidence intervals as shaded areas
#'    \code{ci.plot = "polygon"} or as lines \code{ci.plot = "lines"}.
#' @param ci.args Additional arguments to be passed to the function used
#'    to draw confidence interval. Either \code{\link[graphics]{polygon}} or
#'    \code{\link[graphics]{lines}}.
#' @param ... Additional graphical parameters for the drawn function. See
#'    \code{\link[graphics]{par}}.
#'    yshift and yscale can take intercept and beta if set to NULL
#'
#' @export
plot.cgaim <- function(x, select = NULL, ci = NULL,
  ci.plot = c("polygon", "lines"), ci.args = list(), add = FALSE,
  xcenter = FALSE, xscale = FALSE, yshift = FALSE, yscale = FALSE, ...)
{
  p <- ncol(x$gfit)
  if (is.null(select)){
    select <- seq_len(p)
  } else {
    if (is.character(select)){
      selmatch <- match(select, colnames(x$gfit))
      nas <- is.na(selmatch)
      if (any(nas)) warning(paste0("Incorrect names removed: ",
        paste(select[nas], collapse = ", ")))
      select <- selmatch[!nas]
    } else {
      inc <- select > p
      if (any(inc)) warning(paste0("Incorrect indices removed: ",
        paste(select[inc], collapse = ", ")))
      select <- select[!inc]
    }
  }
  nsel <- length(select)
  if (nsel > 1) {
    if (add) {
      warning("'add = T' should be used with 'select' to add a single smooth")
    } else {
      grDevices::devAskNewPage(TRUE) 
    }
  }
  plotfun <- ifelse(add, graphics::lines, graphics::plot)
  # Initialize centering and scaling
  if (is.numeric(xcenter)) xcenter <- rep_len(xcenter, nsel)
  if (is.numeric(xscale)) xscale <- rep_len(xscale, nsel)
  if(isTRUE(yshift)) yshift <- x$beta[1]
  if(isFALSE(yshift)) yshift <- 0
  yshift <- rep_len(yshift, nsel)
  if(isTRUE(yscale)) yscale <- x$beta[select + 1]
  if(isFALSE(yscale)) yscale <- 1
  yscale <- rep_len(yscale, nsel)
  # Initialize x
  xs <- cbind(x$indexfit, x$covariates)
  xs <- scale(xs[, select], center = xcenter, scale = xscale)
  # Initialize y
  n <- nrow(x$gfit)
  ys <- x$gfit[, select, drop = F] * 
    matrix(yscale, nrow = n, ncol = nsel, byrow = T) + 
    matrix(yshift, nrow = n, ncol = nsel, byrow = T)
  # Initialize parameters
  defParams <- list(ylab = "g", type = "l")
  dots <- list(...)
  # Initialize cis
  if (!is.null(ci)){
    ci.plot <- match.arg(ci.plot)
    allcis <- ci$g[,select,,drop = F]
    allcis[,,1] <- allcis[,,1] * 
      matrix(yscale, nrow = n, ncol = nsel, byrow = T) + 
      matrix(yshift, nrow = n, ncol = nsel, byrow = T)
    allcis[,,2] <- allcis[,,2] * 
      matrix(yscale, nrow = n, ncol = nsel, byrow = T) + 
      matrix(yshift, nrow = n, ncol = nsel, byrow = T)
    if (ci.plot == "polygon"){
      defArgs <- list(border = NA, col = "grey")      
    }
    if (ci.plot == "lines"){
      defArgs <- list(lty = 2)
    }
    ci.args <- modifyList(ci.args, defArgs)
  }
  # Loop to plot
  for (j in seq_len(nsel)){
    xy <- cbind(xs[,j], ys[,j])
    jord <- order(xy[,1])
    xy <- xy[jord,]
    defParams$xlab <- colnames(xs)[j]
    jpars <- modifyList(dots, defParams)
    jpars$x <- xy[,1]
    jpars$y <- xy[,2]
    if (!is.null(ci)){
      cixy <- cbind(xy[,1], allcis[jord,j,])
      if (is.null(jpars$ylim)){
        jpars$ylim <- range(allcis[,j,])
      }
      do.call(plotfun, jpars)
      if (ci.plot == "polygon"){
        ci.args$x <- c(cixy[,1], rev(cixy[,1]))
        ci.args$y <- c(cixy[,2], rev(cixy[,3]))
        do.call(graphics::polygon, ci.args)
        do.call(graphics::points, jpars)
      }
      if (ci.plot == "lines"){
        ci.args$x <- cixy[,1]
        ci.args$y <- cixy[,2]
        do.call(graphics::lines, ci.args)
        ci.args$y <- cixy[,3]
        do.call(graphics::lines, ci.args)
      }
    } else {
      do.call(plotfun, jpars)
    }
  }
  grDevices::devAskNewPage(FALSE)
  invisible()  
}