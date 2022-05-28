vcov_beta <- function(object){
  # Compute unscaled covariance matrix
  X <- cbind(1, object$gfit)
  unsc.vcov <- qr.solve(crossprod(X * sqrt(object$weights)))
  
  # Sigma estimator
  sigsq <- residual_sd(object)
  
  # Compute vcov matrix
  unsc.vcov * sigsq
}
