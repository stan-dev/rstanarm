# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code

set.seed(12345)

MODELS_HOME <- "exec"
fsep <- .Platform$file.sep
if (!file.exists(MODELS_HOME)) {
  MODELS_HOME <- sub(paste0("tests.*", fsep, "testthat$"), 
                     paste0("rstanarm", fsep, "exec"), getwd())
}
if (!file.exists(MODELS_HOME)) {
  MODELS_HOME <- sub(paste0("tests.*", fsep, "testthat$"), "exec", getwd())
}
if (!file.exists(MODELS_HOME)) {
  MODELS_HOME <- system.file("exec", package = "rstanarm") 
}

context("setup")
test_that("Stan programs are available", {
  message(MODELS_HOME)
  expect_true(file.exists(MODELS_HOME))
  expect_true(file.exists(file.path(system.file("chunks", package = "rstanarm"), 
                                    "common_functions.stan")))
  
})
  
library(rstan)
Sys.unsetenv("R_TESTS")

functions <- sapply(dir(MODELS_HOME, pattern = "stan$", full.names = TRUE), function(f) {
  # mc <- scan(file = f, what = "character", sep = "\n", quiet = TRUE)
  mc <- readLines(f)
  start <- grep("^functions[[:blank:]]*\\{[[:blank:]]*$", mc)
  if (length(start) == 1) {
    end <- grep("^}[[:blank:]]*$", mc)[1]
    return(mc[(start + 1L):(end - 1L)])
  }
  else return(as.character(NULL))
})
print(MODELS_HOME)
functions <- c(readLines(file.path(system.file("chunks", package = "rstanarm"), 
                                   "common_functions.stan")), unlist(functions))
model_code <- paste(c("functions {", functions, "}", "model {}"), collapse = "\n")
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

# Negative Binomial
context("Negative Binomial")
test_that("pw_nb return expected results", {
  for (i in 1:length(links)) {
    y <- sample.int(10, size = N, replace = TRUE)
    eta <- abs(rnorm(N))
    linkinv <- poisson(link = links[i])$linkinv
    theta <- rexp(1)
    ll <- dnbinom(y, mu = linkinv(eta), size = theta, log = TRUE)
    expect_true(all.equal(ll,  pw_nb(y, eta, theta, i)), info = links[i])
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

# Gamma GLM
test_that("linkinv_gamma returns expected results", {
  for (i in 1:length(links)) {
    eta <- rexp(N)
    linkinv <- Gamma(link = links[i])$linkinv
    expect_true(all.equal(linkinv(eta), linkinv_gamma(eta, i)), info = links[i])
  }
})
test_that("pw_gamma returns expected results", {
  for (i in 1:length(links)) {
    eta <- rexp(N)
    shape <- rexp(1)
    linkinv <- Gamma(link = links[i])$linkinv
    y <- rgamma(N, shape, rate = 1 / linkinv(eta))
    expect_true(all.equal(dgamma(y, shape = shape, rate = shape / linkinv(eta), log = TRUE),
                          pw_gamma(y, eta, shape, i)), info = links[i])
  }
})
test_that("pw_gamma implies an actual density", {
  for (i in 1:length(links)) {
    eta <- rexp(1)
    shape <- rexp(1)
    foo <- function(y) {
      exp(pw_gamma(y, rep(eta, length(y)), shape, i))
    }
    expect_true(all.equal(1, integrate(foo, lower = 0, upper = Inf)$value, tol = 1e-5))    
  }
})
test_that("GammaReg_log returns the expected results", {
  for (i in 1:length(links)) {
    eta <- rexp(N)
    shape <- rexp(1)
    linkinv <- Gamma(link = links[i])$linkinv
    y <- rgamma(N, shape, rate = 1 / linkinv(eta))
    expect_true(all.equal(sum(dgamma(y, shape = shape, 
                                     rate = shape / linkinv(eta), log = TRUE)),
                          GammaReg_log(y, eta, shape, i, sum(log(y)))), info = links[i])
  }
})
  
# Inverse Gaussian GLM
links <- c(links, "1/mu^2")
test_that("linkinv_inv_gaussian returns expected results", {
  for (i in 1:length(links)) {
    eta <- rgamma(N, 2, 1)
    linkinv <- inverse.gaussian(link = links[i])$linkinv
    expect_true(all.equal(linkinv(eta), linkinv_inv_gaussian(eta, i)), info = links[i])
  }
})
rinvGauss <- function(n, mu, lambda) {
  # from https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
  y <- rchisq(n, 1)
  mu2 <- mu^2
  x <- mu + 0.5 * mu2 * y / lambda - 0.5 * mu / lambda *
       sqrt(4 * mu * lambda * y + mu2 * y^2)
  z <- runif(n)
  out <- ifelse(z <= mu / (mu + x), x, mu2 / x)
  return(out)
}
dinvGauss <- function(x, mu, lambda, log = FALSE) {
  out <- 0.5 * log(0.5 * lambda / pi) - 1.5 * log(x) - 
         0.5 * lambda / mu^2 * (x - mu)^2 / x
  if (!log) out <- exp(out)
  return(out)
}
test_that("pw_inv_gaussian returns expected results", {
  for (i in 1:length(links)) {
    eta <- rgamma(N, 2, 1)
    lambda <- rexp(1)
    linkinv <- inverse.gaussian(link = links[i])$linkinv
    y <- rinvGauss(N, linkinv(eta), lambda)
    expect_true(all.equal(dinvGauss(y, linkinv(eta), lambda, log = TRUE),
                          pw_inv_gaussian(y, eta, lambda, i, log(y), sqrt(y))), info = links[i])
  }
})
test_that("pw_inv_gaussian implies an actual density", {
  for (i in 1:length(links)) {
    eta <- rgamma(1, 2, 1)
    lambda <- rexp(1)
    foo <- function(y) {
      exp(pw_inv_gaussian(y, rep(eta, length(y)), lambda, i, log(y), sqrt(y)))
    }
    expect_true(all.equal(1, integrate(foo, lower = 0, upper = Inf)$value, tol = 1e-4))    
  }
})
test_that("inv_gaussian returns expected results", {
  for (i in 1:length(links)) {
    eta <- rgamma(N, 2, 1)
    lambda <- rexp(1)
    linkinv <- inverse.gaussian(link = links[i])$linkinv
    y <- rinvGauss(N, linkinv(eta), lambda)
    expect_true(all.equal(sum(dinvGauss(y, linkinv(eta), lambda, log = TRUE)),
                          inv_gaussian_log(y, linkinv_inv_gaussian(eta,i), 
                                           lambda, sum(log(y)), sqrt(y))), 
                info = links[i])
  }
})

# lm
N <- 99L
context("lm")
test_that("ll_mvn_ols_qr_lp returns expected results", {
  X <- matrix(rnorm(2 * N), N, 2)
  y <- 1 + X %*% c(2:3) + rnorm(N)
  ols <- lm.fit(cbind(1,X), y)
  b <- coef(ols)
  X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = "-")
  intercept <- 0.5
  beta <- rnorm(2)
  sigma <- rexp(1)
  SSR <- crossprod(residuals(ols))[1]
  ll <- sum(dnorm(y, intercept + X %*% beta, sigma, log = TRUE))
  decomposition <- qr(X)
  Q <- qr.Q(decomposition)
  R <- qr.R(decomposition)
  R_inv <- qr.solve(decomposition, Q)
  b <- R %*% b[-1]
  beta <- R %*% beta
  expect_true(all.equal(ll, ll_mvn_ols_qr_lp(beta, b, intercept, mean(y), 
                                             SSR, sigma, N)))
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
    expect_equal(log(Pr), pw_polr(y, eta, zeta, i, 1), info = links[i])
  }
})
rdirichlet <- function(n, alpha) {
  # from MCMCpack::rdirichlet and licensed under the GPL
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}
context("polr")
test_that("make_cutpoints returns expected results", {
  J <- 5L
  for (i in 1:length(links)) {
    p <- rdirichlet(1, rep(1,J))[1,]
    cutpoints <- make_cutpoints(p, 1, i)
    for (j in 1:length(cutpoints)) {
      expect_true(all.equal(sum(p[1:j]), CDF_polr(cutpoints[j], i)))
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

context("glmer")
test_that("the Stan equivalent of lme4's Z %*% b works", {
  stopifnot(require(lme4))
  stopifnot(require(Matrix))
  test_lme4 <- function(group) {
    Lambdati <- group$Lambdat
    Lind <- group$Lind
    theta <- group$theta
    
    group <- rstanarm:::pad_reTrms(Z = t(group$Zt), cnms = group$cnms, 
                                   flist = group$flist)
    Z <- group$Z
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attr(group$flist, "assign"), function(i) nlevels(group$flist[,i]))
    
    len_theta_L <- sum(choose(p,2), p)
    expect_true(len_theta_L == length(theta))
    dispersion <- runif(1)
    tau <- as.array(rgamma(length(p), shape = 1, scale = 1))
    scale <- as.array(abs(rcauchy(length(p))))
    zeta <- as.array(rgamma(sum(p[p > 1]), shape = 1, scale = 1))
    rho <- as.array(rbeta(sum(p - 1), 1, 1))
    z_T <- as.array(rnorm(sum(pmax(0, choose(p,2) - 1))))
    
    theta_L <- make_theta_L(len_theta_L, p, dispersion, tau, scale, zeta, rho, z_T)
    expect_true(all(theta_L[theta == 1] > 0))
    Lambdati@x <- theta_L[Lind]
    
    z_b <- rnorm(ncol(Z))
    b <- make_b(z_b, theta_L, p, l)
    mark <- colnames(Z) == ""
    expect_equal(b[!mark], as.vector(Matrix::t(Lambdati) %*% z_b[!mark]), 
                 tol = 1e-14)
    
    parts <- extract_sparse_parts(Z)
    Zb <- test_csr_matrix_times_vector(nrow(Z), ncol(Z), parts$w, 
                                       parts$v, parts$u, b)
    expect_equal(Zb, as.vector(Z %*% b), tol = 1e-14)
  }
  test_lme4(glFormula(Reaction ~ Days + (Days | Subject), data = sleepstudy)$reTrms)
  test_lme4(glFormula(Reaction ~ Days + (Days || Subject), data = sleepstudy)$reTrms)
  test_lme4(glFormula(Reaction ~ Days + (1 | Subject), data = sleepstudy)$reTrms)
  test_lme4(glFormula(cbind(incidence, size - incidence) ~ period + (1 | herd),
                            data = cbpp, family = binomial)$reTrms)
  cbpp$obs <- 1:nrow(cbpp)
  test_lme4(glFormula(cbind(incidence, size - incidence) ~ period +
                        (1 | herd) +  (1|obs), family = binomial, data = cbpp)$reTrms)
  data(toenail, package = "HSAUR3")
  test_lme4(glFormula(outcome ~ visit + treatment + (visit|treatment) + (1|patientID),
                      data=toenail, family = binomial)$reTrms)
  data(clouds, package = "HSAUR3")
  test_lme4(glFormula(rainfall ~ sne + cloudcover + prewetness + echomotion + 
                        (1 + sne + cloudcover + prewetness|seeding) +  
                        (1 + sne + cloudcover + prewetness||echomotion),
                      data=clouds, family = gaussian)$reTrms)
  test_lme4(glFormula(angle ~ recipe + temp + (1|recipe:replicate), data = cake)$reTrms)
})

