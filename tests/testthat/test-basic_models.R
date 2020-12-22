#----------------------------
# Create data
#----------------------------

#---- A very simple single-index model
set.seed(1989)
n <- 200
x1 <- rnorm(n)
x2 <- x1 + rnorm(n)
z <- x1 + x2
y <- z + rnorm(n)
df1 <- data.frame(y, x1, x2)

#---- A simple two-index model
set.seed(2020)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rnorm(n)
mu <- 4 * exp(8 * x1) / (1 + exp(8 * x1)) + exp(x3)
y <- mu + rnorm(n)
df2 <- data.frame(y, x1, x2, x3, x4)

#----------------------------
# Test models without constraints
#----------------------------
ans <- cgaim(y ~ g(x1, x2), data = df1) 
ans <- cgaim(y ~ g(x1, x2) + g(x3, x4), data = df2)

#----------------------------
# Test alpha constraints
#----------------------------

#----- Monotonicity constraints
# Single-index
ans1 <- cgaim(y ~ g(x1, x2, acons = list(monotone = 1)), 
  data = df1)
ans2 <- cgaim(y ~ g(x1, x2, acons = list(monotone = -1)), 
  data = df1)
# Two-index
ans3 <- cgaim(y ~ g(x1, x2, acons = list(monotone = 1)) + 
    g(x3, x4, acons = list(monotone = 1)), 
  data = df2)
ans4 <- cgaim(y ~ g(x1, x2, acons = list(monotone = -1)) + 
    g(x3, x4, acons = list(monotone = -1)), 
  data = df2)
  
test_that("monotonicity constraints on alpha work", {
  expect_gte(ans1$alpha.control$Cmat %*% unlist(ans1$alpha), 0)
  expect_gte(ans2$alpha.control$Cmat %*% unlist(ans2$alpha), 0)
  expect_true(all(ans3$alpha.control$Cmat %*% unlist(ans3$alpha) > 0))
  expect_true(all(ans4$alpha.control$Cmat %*% unlist(ans4$alpha) > 0))
})

#----- Sign constraints
# Single-index
ans1 <- cgaim(y ~ g(x1, x2, acons = list(sign.const = 1)), 
  data = df1)
ans2 <- cgaim(y ~ g(x1, x2, acons = list(sign.const = -1)), 
  data = df1)
# Two-index
ans3 <- cgaim(y ~ g(x1, x2, acons = list(sign.const = 1)) + 
    g(x3, x4, acons = list(sign.const = 1)), 
  data = df2)
ans4 <- cgaim(y ~ g(x1, x2, acons = list(sign.const = -1)) + 
    g(x3, x4, acons = list(sign.const = -1)), 
  data = df2)
  
test_that("sign constraints on alpha work", {
  expect_true(all(ans1$alpha.control$Cmat %*% unlist(ans1$alpha) > 0))
  expect_true(all(ans2$alpha.control$Cmat %*% unlist(ans2$alpha) > 0))
  expect_true(all(ans3$alpha.control$Cmat %*% unlist(ans3$alpha) > 0))
  expect_true(all(ans4$alpha.control$Cmat %*% unlist(ans4$alpha) > 0))
})

#----------------------------
# Test smoothing constraints
#----------------------------

#----- Monotone increasing constraints
ans1 <- cgaim(y ~ g(x1, x2, fcons = "inc"), 
  data = df1, smooth_method = "scam")
ans2 <- cgaim(y ~ g(x1, x2, fcons = "inc") + g(x3, x4, fcons = "inc"), 
  data = df2, smooth_method = "scam")

test_that("Monotone increasing constraint on smooths works with 'scam'",{
  expect_true(all(diff(ans1$gfit[order(ans1$indexfit)]) >= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,1]), 1]) >= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,2]), 2]) >= 0))
})

ans1 <- cgaim(y ~ g(x1, x2, fcons = "inc"), 
  data = df1, smooth_method = "scar")
ans2 <- cgaim(y ~ g(x1, x2, fcons = "inc") + g(x3, x4, fcons = "inc"), 
  data = df2, smooth_method = "scar")

test_that("Monotone increasing constraint on smooths works with 'scar'",{
  expect_true(all(diff(ans1$gfit[order(ans1$indexfit)]) >= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,1]), 1]) >= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,2]), 2]) >= 0))
})

# ans1 <- cgaim(y ~ g(x1, x2, fcons = "inc"), 
#   data = df1, smooth_method = "cgam")
# ans2 <- cgaim(y ~ g(x1, x2, fcons = "inc") + g(x3, x4, fcons = "inc"), 
#   data = df2, smooth_method = "cgam")
# 
# test_that("Monotone increasing constraint on smooths works with 'cgam'",{
#   expect_true(all(diff(ans1$gfit[order(ans1$indexfit)]) >= 0))
#   expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,1]), 1]) >= 0))
#   expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,2]), 2]) >= 0))
# })

#----- Monotone decreasing constraints
ans1 <- cgaim(y ~ g(x1, x2, fcons = "dec"), 
  data = df1, smooth_method = "scam")
ans2 <- cgaim(y ~ g(x1, x2, fcons = "dec") + g(x3, x4, fcons = "dec"), 
  data = df2, smooth_method = "scam")

test_that("Monotone increasing constraint on smooths works with 'scam'",{
  expect_true(all(diff(ans1$gfit[order(ans1$indexfit)]) <= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,1]), 1]) <= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,2]), 2]) <= 0))
})

ans1 <- cgaim(y ~ g(x1, x2, fcons = "dec"), 
  data = df1, smooth_method = "scar")
ans2 <- cgaim(y ~ g(x1, x2, fcons = "dec") + g(x3, x4, fcons = "dec"), 
  data = df2, smooth_method = "scar")

test_that("Monotone increasing constraint on smooths works with 'scar'",{
  expect_true(all(diff(ans1$gfit[order(ans1$indexfit)]) <= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,1]), 1]) <= 0))
  expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,2]), 2]) <= 0))
})

# ans1 <- cgaim(y ~ g(x1, x2, fcons = "dec"), 
#   data = df1, smooth_method = "cgam")
# ans2 <- cgaim(y ~ g(x1, x2, fcons = "dec") + g(x3, x4, fcons = "dec"), 
#   data = df2, smooth_method = "cgam")
# 
# test_that("Monotone increasing constraint on smooths works with 'cgam'",{
#   expect_true(all(diff(ans1$gfit[order(ans1$indexfit)]) <= 0))
#   expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,1]), 1]) <= 0))
#   expect_true(all(diff(ans2$gfit[order(ans2$indexfit[,2]), 2]) <= 0))
# })
