#' Predictions using the GAIM model
#'
#' @param object A \code{gaim} object resulting from a call to \code{gaim}.
#' @param newdata A list containing the new data to predict. Must be a list
#'    similar to the one passed in \code{x} in the call to \code{gaim}.
#'    If missing, uses the training data.  
predict.gaim <- function(object, newdata){
  #! Add checks that it is indeed a gaim object
  if (missing(newdata) || is.null(newdata)){
    x <- object$x
  } else {
    x <- newdata
  }
  x <- lapply(x, as.matrix) # For the matrix product
  alpha <- object$alpha
  z <- mapply("%*%", x, alpha)
  z <- as.data.frame(z)
  varnames <- attr(object$am.fit$terms, "term.labels")
  if (!all.equal(colnames(z), varnames)){
    warning("Variable names do not correspond")
    colnames(z) <- varnames
  } 
  #! Don't forget to chack for the indices magnitude
  linear.predictor <- predict(object$am.fit, newdata = z)
  yhat <- linear.predictor + object$coef[1]
  return(yhat)
}