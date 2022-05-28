#----- Simple model
set.seed(2020)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rnorm(n)
x5 <- rnorm(n)
alpha <- list(c(.8, .2), c(.5, .5))
mu <- 4 * 
  exp(cbind(x1, x2) %*% alpha[[1]]) / (1 + exp(cbind(x1, x2) %*% alpha[[1]])) + 
  exp(x3) + sin(cbind(x4, x5) %*% alpha[[2]])

y <- mu + rnorm(n)
df3 <- data.frame(y, x1, x2, x3, x4, x5)

ans <- cgaim(y ~ g(x1, x2, label = "i1", s_opts = list(sp = 0)) + 
    s(x3, s_opts = list(bs = "cr")) + g(x4, x5, label = "i2", fcons = "inc"),
  data = df3)

#----- Evaluate coverage

ns <- 100
B <- 1000

# Simulate
ysim <- drop(mu) + matrix(rnorm(n * ns), n, ns)

# Loop
cl <- makeCluster(max(1, detectCores() - 2))
registerDoParallel(cl)
res <- foreach(i = seq_len(ns), .packages = "cgaim", .combine = rbind) %dopar% {
  dfi <- data.frame(y = ysim[,i], x1, x2, x3, x4, x5)
  ans <- cgaim(y ~ g(x1, x2, label = "i1", s_opts = list(sp = 0)) + 
      s(x3, s_opts = list(sp = 0)) + 
      g(x4, x5, label = "i2", s_opts = list(sp = 0)),
    data = dfi)
  ci <- confint(ans, parm = "alpha", B = B, type = "norm")
  data.frame(simu = i, alpha = 1:4, est = unlist(ans$alpha), ci$alpha,
    sd = diag(cgaim:::vcov_alpha(ans)))
}
stopCluster(cl)

# Get coverage
covered <- data.table::between(rep(unlist(alpha), ns), res[,4], res[,5])
aggregate(covered ~ alpha, res, mean)

# Get coverage 2
covered <- data.table::between(rep(unlist(alpha), ns), 
  res[,3] - 1.96 * sqrt(res[,6]), res[,3] + 1.96 * sqrt(res[,6]))
aggregate(covered ~ alpha, res, mean)

# Get vcov
estmat <- matrix(res$est, nrow = ns, ncol = 4, byrow = T)
estvcov <- var(estmat)


#----- Test with nls
x <- 1:10
y <- 2*x + 3   
set.seed(27)
yeps <- y + rnorm(length(y), sd = 0.01) # added noise
res <- nls(yeps ~ a + b*x, start = list(a = 0.12345, b = 0.54321))

nlsres <- replicate(ns, {
  yeps <- y + rnorm(length(y), sd = 2) # added noise
  res <- nls(yeps ~ a + b*x)
  v <- diag(vcov(res))
  cbind(coef(res) - 1.96 * sqrt(v), coef(res) + 1.96 * sqrt(v))
})

# Coverage
covered <- apply(nlsres, 3, function(x) data.table::between(c(3, 2), x[,1], x[,2]))
rowMeans(covered)
