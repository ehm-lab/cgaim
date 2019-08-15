library(MASS)

X <- replicate(2, mvrnorm(100, rep(0, 2), diag(2)), simplify = F)
Alpha <- list(c(.5, .5), c(.9, .1))
Z <- mapply("%*%", X, Alpha)

test_that("predict training is identical to fitted values", {
  # Simple case without noise
  Y <- exp(Z[,1]) + sin(Z[,2])
  res_gaim <- aim(Y, X)
  
  expect_equal(predict.gaim(res_gaim, newdata = X), res_gaim$fitted,
    check.attributes = FALSE)
    
  # With noise
  Y <- exp(Z[,1]) + sin(Z[,2]) + rnorm(100, 0, .2)
  res_gaim <- aim(Y, X)
  
  expect_equal(predict.gaim(res_gaim, newdata = X), res_gaim$fitted,
    check.attributes = FALSE)
  
  # With different magnitudes
  Y <- 5 * exp(Z[,1]) + .5 * sin(Z[,2]) + rnorm(100, 0, .2)
  res_gaim <- aim(Y, X)
  
  expect_equal(predict.gaim(res_gaim, newdata = X), res_gaim$fitted,
    check.attributes = FALSE)
})


test_that("predict validation data work", {
  Y <- exp(Z[,1]) + sin(Z[,2]) + rnorm(100, 0, .2)
  res_gaim <- aim(Y, X)
  newx <- replicate(2, mvrnorm(5, rep(0,2), diag(2)), simplify = F)
  yhat <- predict.gaim(res_gaim, newdata = newx)
  
  expect_length(yhat, 5)
})