#' Predictions from a fitted CGAIM object
#'
#' Uses a fitted \code{cgaim} object and computes prediction for the
#'    observed data or new data. Predicts the response, indices or 
#'    ridge functions values at the provided data.
#' 
#' When \code{newdata} is provided, it must contain all variables used in
#'    the fitted model.
#'
#' @param object A \code{gaim} object.
#' @param newdata A list or data.frame containing the new data to predict.
#'    If missing, fitted values from the model are returned.
#' @param type A character indicating the type of prediction to return.
#'    When \code{type = "response"}, the predicted response is returned. When
#'    \code{type = "terms"}, returns each ridge and smooth function separately.
#'    When \code{type = "indices"}, returns predicted indices values.
#' @param select A numeric or character vector indicating terms to return
#'    when \code{type = "terms"} or \code{type = "indices"}.
#' @param na.action A function indicating how to treat NA values when 
#'    when \code{newdata} is not missing.
#' @param ... For compatibility with the default \code{predict} method. Unused
#'    at the moment.
#'
#' @return When \code{type = "response"} returns a vector of predicted response.
#'    When \code{type = "terms"} returns a d-column matrix where d is the 
#'    number of ridge and smooth terms. When \code{type = "indices"}, returns
#'    a p-column matrix where p is the number of created indices.
#'
#' @export
predict.cgaim <- function(object, newdata, 
  type = c("response", "terms", "indices"), select = NULL, 
  na.action = "na.pass", ...)
{
  type <- match.arg(type)
  n <- length(object$fitted)
  p <- length(object$beta) - 1
  mt <- stats::terms(object)
  if (!missing(newdata)){
    mt <- stats::delete.response(mt)
    gind <- attr(mt, "specials")$g
    mfind <- stats::model.frame(mt[gind], data = newdata, 
      na.action = na.action)
    alphas <- object$alpha
    newindex <- mapply("%*%", mfind, alphas)
    colnames(newindex) <- names(alphas)
    if (type == "indices"){
      if (!is.null(select)) newindex <- newindex[,select, drop = F]
      return(newindex)
    }
    Xterms <- c(data.frame(newindex), newdata[names(object$covariates)])
    sind <- attr(mt, "specials")$s
    objx <- cbind(object$indexfit, object$covariates)
    gterms <- matrix(0, nrow(mfind), p)
    for (j in 1:p){
      if (is.numeric(Xterms[[j]])){
        inter_fun <- ifelse(j %in% c(gind, sind), "spline", "approx")
        gterms[,j] <- suppressWarnings(do.call(inter_fun, 
          list(x = objx[,j], y = object$gfit[,j], xout = Xterms[[j]]))$y) 
      } else {
        faccoefs <- unique(object$gfit[,j])
        names(faccoefs) <- unique(objx[,j])
        gterms[,j] <- faccoefs[Xterms[[j]]]
      }
    }
    colnames(gterms) <- colnames(object$gfit)
    if (type == "terms"){
      if (!is.null(select)) gterms <- gterms[,select, drop = F]
      return(gterms)
    }
    betas <- object$beta
    yhat <- cbind(1, gterms) %*% betas
    return(yhat)
  } else {
    out <- switch(type,
      response = object$fitted,
      terms = object$gfit * matrix(object$beta[-1], n, p, byrow = T),
      indices = object$indexfit
    )
    return(out)
  }
}
