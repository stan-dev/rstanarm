# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code

MODELS_HOME <- "exec"
fsep <- .Platform$file.sep
if (!dir.exists(MODELS_HOME)) {
  MODELS_HOME <- sub(paste0("tests", fsep, "testthat$"), 
                     paste0("rstanarm", fsep, "exec"), getwd())
}
if (!dir.exists(MODELS_HOME)) {
  MODELS_HOME <- sub(paste0("tests", fsep, "testthat$"), "exec", getwd())
}

context("setup")
test_that("Stan programs are available", {
  message(MODELS_HOME)
  expect_true(dir.exists(MODELS_HOME))  
})
  
stopifnot(require(rstan))
Sys.unsetenv("R_TESTS")

functions <- sapply(dir(MODELS_HOME, pattern = "stan$", full.names = TRUE), function(f) {
  mc <- scan(text = stanc(f)$model_code, what = "character", sep = "\n",
             quiet = TRUE)
  start <- grep("^functions[[:blank:]]*\\{[[:blank:]]*$", mc)
  if (length(start) == 1) {
    end <- grep("^}[[:blank:]]*$", mc)[1]
    return(mc[(start + 1L):(end - 1L)])
  }
  else return(as.character(NULL))
})

model_code <- paste(c("functions {", unlist(functions), "}", "model {}"),
                    collapse = "\n")
expose_stan_functions(stanc(model_code = model_code, model_name = "Stan Functions"))
N <- 99L

# bernoulli
links <- c("logit", "probit", "cauchit", "log", "cloglog")

context("Bernoulli")
test_that("linkinv_bern returns expected results", {
  for (i in 1:length(links)) {
    eta <- -abs(rnorm(N))
    linkinv <- binomial(link = links[i])$linkinv
    expect_true(all.equal(linkinv(eta), 
                          linkinv_bern(eta, i)), info = links[i])
  }
})
context("Bernoulli")
test_that("pw_bern and ll_bern_lp return expected results", {
  for (i in 1:length(links)) {
    eta0 <- -abs(rnorm(N))
    eta1 <- -abs(rnorm(N))
    linkinv <- binomial(link = links[i])$linkinv
    ll0 <- dbinom(0, size = 1, prob = linkinv(eta0), log = TRUE)
    expect_true(all.equal(ll0, pw_bern(0, eta0, i)), info = links[i])
    ll1 <- dbinom(1, size = 1, prob = linkinv(eta1), log = TRUE)
    expect_true(all.equal(ll1, pw_bern(1, eta1, i)), info = links[i])
    expect_true(all.equal(sum(ll0, ll1), 
                          ll_bern_lp(eta0, eta1, i, c(N,N))), 
                info = links[i])
  }
})
context("Bernoulli")
test_that("make_upper_bernoulli returns expected results", {
  for (i in 1:length(links)) {
    X0 <- matrix(rnorm(2 * N), N, 2)
    X1 <- matrix(rnorm(2 * N), N, 2)
    beta <- rnorm(2)
    has_offset <- 0L
    offset0 <- matrix(0, nrow = N, ncol = 0)
    offset1 <- matrix(0, nrow = N, ncol = 0)
    if (i == 4) expect_true(is.finite(make_upper_bernoulli(i, 
                            X0, X1, beta, has_offset, offset0, offset1)))
    else expect_true(Inf == make_upper_bernoulli(i, 
                     X0, X1, beta, has_offset, offset0, offset1), info = links[i])
  }
})

context("Bernoulli")
test_that("dp_deta returns expected results", {
  for (i in 1:length(links)) {
    eta <- rnorm(1)
    expect_true(all.equal(make.link(links[i])$mu.eta(eta), dp_deta(eta, i)))
  }
})

# Binomial
trials <- 10L
context("Binomial")
test_that("linkinv_binom returns expected results", {
  for (i in 1:length(links)) {
    eta <- -abs(rnorm(N))
    linkinv <- binomial(link = links[i])$linkinv
    expect_true(all.equal(linkinv(eta), 
                          linkinv_binom(eta, i)), info = links[i])
  }
})
context("Bernoulli")
test_that("pw_binom and ll_binom_lp return expected results", {
  for (i in 1:length(links)) {
    eta <- -abs(rnorm(N))
    y <- sample.int(trials, size = N, replace = TRUE)
    linkinv <- binomial(link = links[i])$linkinv
    ll <- dbinom(y, size = trials, prob = linkinv(eta), log = TRUE)
    expect_true(all.equal(ll,  pw_binom(y, rep(trials, N), eta, i)), info = links[i])
    expect_true(all.equal(sum(ll), ll_binom_lp(y, rep(trials, N), eta, i) + 
                ifelse(i > 3, sum(lchoose(trials, y)), 0)), info = links[i])
  }
})
context("Bernoulli")
test_that("make_upper_bernoulli returns expected results", {
  for (i in 1:length(links)) {
    X <- matrix(rnorm(2 * N), N, 2)
    beta <- rnorm(2)
    has_offset <- 0L
    offset <- matrix(0, nrow = N, ncol = 0)
    if (i == 4) expect_true(is.finite(make_upper_binomial(i, 
                            X, beta, has_offset, offset)))
    else expect_true(Inf == make_upper_binomial(i, 
                            X, beta, has_offset, offset), info = links[i])
  }
})

# Count GLM
links <- c("log", "identity", "sqrt")

context("Poisson")
test_that("linkinv_count returns expected results", {
  for (i in 1:length(links)) {
    eta <- abs(rnorm(N))
    linkinv <- poisson(link = links[i])$linkinv
    expect_true(all.equal(linkinv(eta), 
                          linkinv_count(eta, i)), info = links[i])
  }
})
context("Poisson")
test_that("pw_pois return expected results", {
  for (i in 1:length(links)) {
    y <- sample.int(10, size = N, replace = TRUE)
    eta <- abs(rnorm(N))
    linkinv <- poisson(link = links[i])$linkinv
    ll <- dpois(y, linkinv(eta), log = TRUE)
    expect_true(all.equal(ll,  pw_pois(y, eta, i)), info = links[i])
  }
})
context("Poisson")
test_that("make_lower_count returns expected results", {
  for (i in 1:length(links)) {
    X <- matrix(rnorm(2 * N), N, 2)
    beta <- rnorm(2)
    has_offset <- 0L
    offset <- matrix(0, nrow = N, ncol = 0)
    if (i == 2) expect_true(is.finite(make_lower_count(i, 
                                      X, beta, has_offset, offset)))
    else expect_true(-Inf == make_lower_count(i, 
                             X, beta, has_offset, offset), info = links[i])
  }
})

# Gaussian GLM
links <- c("identity", "log", "inverse")

context("Gaussian")
test_that("linkinv_gauss returns expected results", {
  for (i in 1:length(links)) {
    eta <- rnorm(N)
    linkinv <- gaussian(link = links[i])$linkinv
    expect_true(all.equal(if (i == 2) eta else linkinv(eta), 
                          linkinv_gauss(eta, i)), info = links[i])
  }
})
context("Gaussian")
test_that("pw_gauss returns expected results", {
  for (i in 1:length(links)) {
    eta <- rnorm(N)
    linkinv <- gaussian(link = links[i])$linkinv
    if (i == 2)
      expect_true(all.equal(dnorm(0, mean = eta, log = TRUE),
                            pw_gauss(rep(1,N), eta, 1, i)), info = links[i])
    else 
      expect_true(all.equal(dnorm(0, mean = linkinv(eta), log = TRUE),
                            pw_gauss(rep(0,N), eta, 1, i)), info = links[i])
  }
})


# lm
N <- 99L
context("lm")
test_that("ll_mvn_ols_lp returns expected results", {
  X <- matrix(rnorm(2 * N), N, 2)
  y <- 1 + X %*% c(2:3) + rnorm(N)
  ols <- lm.fit(cbind(1,X), y)
  b <- coef(ols)
  X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = "-")
  XtX <- crossprod(X)
  intercept <- 0.5
  beta <- rnorm(2)
  sigma <- rexp(1)
  expect_true(all.equal(sum(dnorm(y, intercept + X %*% beta, sigma, log = TRUE)),
                        ll_mvn_ols_lp(beta, b[-1], XtX, intercept, mean(y),
                                   crossprod(residuals(ols))[1], sigma, N)))
})

# polr
links <- c("logistic", "probit", "loglog", "cloglog", "cauchit")
context("polr")
test_that("CDF_polr returns expected results", {
  for (i in 1:length(links)) {
    x <- rnorm(1)
    if (i == 1) linkinv <- make.link("logit")$linkinv
    else if (i == 3) linkinv <- rstanarm:::pgumbel
    else linkinv <- make.link(links[i])$linkinv
    expect_true(all.equal(linkinv(x), CDF_polr(x, i)))
  }
})
context("polr")
test_that("pw_polr returns expected results", {
  J <- 3
  for (i in 1:length(links)) {
    x <- matrix(rnorm(N * 2), nrow = N, ncol = 2)
    beta <- rnorm(2)
    zeta <- sort(rnorm(J-1))
    eta <- c(x %*% beta)
    y <- apply(rmultinom(N, 1, prob = rep(1/J, J)) == 1, 2, which)
    model <- MASS::polr(as.factor(y) ~ x, method = links[i], 
                        start = c(beta, zeta), control = list(maxit = 0))
    Pr <- fitted(model)
    Pr <- sapply(1:N, FUN = function(i) Pr[i,y[i]])
    expect_true(all.equal(log(Pr), pw_polr(y, eta, zeta, i)))
  }
})
context("polr")
test_that("inv_Phi returns expected results", {
  x <- rnorm(1)
  expect_true(all.equal(x, inv_Phi(pnorm(x))))
})
context("polr")
test_that("make_cutpoints returns expected results", {
  for (i in 1:length(links)) {
    p <- MCMCpack::rdirichlet(1, rep(1,J))[1,]
    cutpoints <- make_cutpoints(p, 1, i)
    for (j in 1:length(cutpoints)) {
      expect_true(all.equal(sum(p[1:j]), CDF(cutpoints[j], i)))
    }
  }
})
context("polr")
test_that("draw_ystar_rng returns expected results", {
  l <- -0.1
  u <-  0.1
  eta <- 0
  for (i in 1:length(links)) {
    draw <- draw_ystar_rng(l, u, eta, i)
    expect_true(draw > l)
    expect_true(draw < u)
  }
})
