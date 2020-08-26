# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017, 2018, 2019 Trustees of Columbia University
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

Sys.setenv(USE_CXX14 = 1)
set.seed(12345)

MODELS_HOME <- "stan_files"
INCLUDE_DIR <- "include"

context("setup")
test_that("Stan programs are available", {
  expect_true(file.exists(MODELS_HOME))
})

library(rstan)
Sys.unsetenv("R_TESTS")
TBB <- system.file("lib", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
SH  <- system.file(ifelse(.Platform$OS.type == "windows", "libs", "lib"), 
                   .Platform$r_arch, package = "StanHeaders",  mustWork = TRUE)
Sys.setenv(LOCAL_LIBS = paste0("-L", shQuote(TBB), " -tbb -tbbmalloc ",
                               "-L", shQuote(SH) , " -lStanHeaders"))
# Sys.setenv(PKG_LIBS = Sys.getenv("LOCAL_LIBS"))
Eigen <- dir(system.file("include", "stan", "math", "prim",
                         package = "StanHeaders", mustWork = TRUE),
             pattern = "Eigen.hpp$", full.names = TRUE, recursive = TRUE)[1]
Sys.setenv(PKG_CXXFLAGS = paste("-include", shQuote(Eigen)))

functions <- sapply(dir(MODELS_HOME, pattern = "stan$", full.names = TRUE), function(f) {
  mc <- readLines(f)
  mc <- grep("^#include", mc, invert = TRUE, value = TRUE)
  start <- grep("^functions[[:blank:]]*\\{[[:blank:]]*$", mc)
  if (length(start) == 1) {
    end <- grep("^}[[:blank:]]*$", mc)[1]
    if (end == (start + 1L)) return(as.character(NULL))
    return(mc[(start + 1L):(end - 1L)])
  } else return(as.character(NULL))
})
names(functions) <- basename(names(functions))
functions$polr.stan <- grep("csr_matrix_times_vector2", 
                            functions$polr.stan, 
                            value = TRUE, fixed = TRUE, invert = TRUE)
functions <- c(unlist(lapply(file.path(MODELS_HOME, "functions", 
                             c("common_functions.stan",
                               "bernoulli_likelihoods.stan",
                               "binomial_likelihoods.stan",
                               "continuous_likelihoods.stan",
                               "count_likelihoods.stan", 
                               "SSfunctions.stan")), 
                      FUN = readLines)), unlist(functions))
model_code <- paste(c("functions {", functions, "}"), collapse = "\n")
stanc_ret <- stanc(model_code = model_code, model_name = "Stan Functions",
                   allow_undefined = TRUE)
expose_stan_functions(stanc_ret, rebuild = TRUE, verbose = TRUE)
Rcpp::sourceCpp(file.path(INCLUDE_DIR, "tests.cpp"), rebuild = TRUE, verbose = TRUE)
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
    expect_true(all.equal(sum(ll), ll_binom_lp(y, rep(trials, N), eta, i), info = links[i]))
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
    expect_true(all.equal(linkinv(eta), linkinv_gauss(eta, i)), info = links[i])
  }
})
context("Gaussian")
test_that("pw_gauss returns expected results", {
  for (i in 1:length(links)) {
    eta <- rnorm(N)
    linkinv <- gaussian(link = links[i])$linkinv
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
                          GammaReg(y, eta, shape, i, sum(log(y)))), info = links[i])
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
                          pw_inv_gaussian(y, eta, lambda, i, log(y), sqrt(y))), 
                info = links[i])
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
                          inv_gaussian(y, linkinv_inv_gaussian(eta,i), 
                                       lambda, sum(log(y)), sqrt(y))), 
                info = links[i])
  }
})

# lm
N <- 99L
context("lm")
test_that("ll_mvn_ols... returns expected results", {
  X <- matrix(rnorm(2 * N), N, 2)
  X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = "-")
  y <- 1 + X %*% c(2:3) + rnorm(N)
  ols <- lm.fit(cbind(1,X), y)
  b <- coef(ols)
  intercept <- 0.5
  beta <- rnorm(2)
  sigma <- rexp(1)
  SSR <- crossprod(residuals(ols))[1]
  ll <- sum(dnorm(y, intercept + X %*% beta, sigma, log = TRUE))
  expect_true(all.equal(ll, ll_mvn_ols(c(intercept, beta), b, 
                                       crossprod(cbind(1, X)), SSR, 
                                       sigma, N)))
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
    log_pr <- pw_polr(y, eta, zeta, i, 1)
    log_Pr <- log(Pr)
    good <- is.finite(log_pr) & is.finite(log_Pr) & log_Pr > -30
    expect_equal(log_Pr[good], log_pr[good], info = links[i], 
                 tolerance = 1e-6)
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

# glmer
context("glmer")
if (require(lme4) && require(HSAUR3)) test_that("the Stan equivalent of lme4's Z %*% b works", {
  stopifnot(require(Matrix))
  test_lme4 <- function(group) {
    Lambdati <- group$Lambdat
    Lind <- group$Lind
    theta <- group$theta
    
    group <- rstanarm:::pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                                   flist = group$flist)
    Z <- group$Z
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attr(group$flist, "assign"), function(i) nlevels(group$flist[[i]]))
    
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
    Zb <- csr_matrix_times_vector2_test(nrow(Z), ncol(Z), parts$w, 
                                        parts$v - 1L, parts$u - 1L, b)
    expect_equal(Zb, as.vector(Z %*% b), tol = 1e-14)
    if (all(sapply(group$cnms, FUN = function(x) {
        length(x) == 1 && x == "(Intercept)"
      })) ) {
      V <- matrix(parts$v, nrow = sum(p), ncol = nrow(Z))
      expect_true(all(V == 
                      t(as.matrix(as.data.frame(make_V(nrow(Z), nrow(V), parts$v - 1L))))))
      expect_equal(Zb, apply(V, 2, FUN = function(v) sum(b[v])))
    }
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
  test_lme4(glFormula(diameter ~ (1|plate) + (1|sample), data = Penicillin)$reTrms)
})

context("glmer")
test_that("the Cornish-Fisher expansion from standard normal to Student t works", {
  df <- exp(1) / pi
  approx_t <- sapply(rnorm(1000), FUN = CFt, df = df)
  expect_true(ks.test(approx_t, "pt", df = df, exact = TRUE)$p.value > 0.05)
})

context("nlmer")
test_that("SSasymp works", {
  Lob.329 <- Loblolly[ Loblolly$Seed == "329", ]
  Asym <- 100
  resp0 <- -8.5
  lrc <- -3.2  
  Phi <- cbind(Asym, resp0, lrc)
  expect_true(all.equal(SSasymp( Lob.329$age, Asym, resp0, lrc ),
                        SS_asymp( Lob.329$age, Phi ), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(Lob.329), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSasymp( Lob.329$age, Asym, resp0, lrc ),
                        SS_asymp( Lob.329$age, Phi ), check.attributes = FALSE))
})

context("nlmer")
test_that("SSasympOff works", {
  CO2.Qn1 <- CO2[CO2$Plant == "Qn1", ]
  Asym <- 32; lrc <- -4; c0 <- 43
  Phi <- cbind(Asym, lrc, c0)
  expect_true(all.equal(SSasympOff(CO2.Qn1$conc, Asym, lrc, c0),
                        SS_asympOff(CO2.Qn1$conc, Phi), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(CO2.Qn1), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSasympOff(CO2.Qn1$conc, Asym, lrc, c0),
                        SS_asympOff(CO2.Qn1$conc, Phi), check.attributes = FALSE))
})

context("nlmer")
test_that("SSasympOrig works", {
  Lob.329 <- Loblolly[ Loblolly$Seed == "329", ]
  Asym <- 100; lrc <- -3.2
  Phi <- cbind(Asym, lrc)
  expect_true(all.equal(SSasympOrig(Lob.329$age, Asym, lrc),
                        SS_asympOrig(Lob.329$age, Phi), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(Lob.329), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSasympOrig(Lob.329$age, Asym, lrc),
                        SS_asympOrig(Lob.329$age, Phi), check.attributes = FALSE))
})

context("nlmer")
test_that("SSbiexp works", {
  Indo.1 <- Indometh[Indometh$Subject == 1, ]
  A1 <- 3; lrc1 <- 1; A2 <- 0.6; lrc2 <- -1.3
  Phi <- cbind(A1, lrc1, A2, lrc2)
  expect_true(all.equal(SSbiexp( Indo.1$time, A1, lrc1, A2, lrc2 ),
                        SS_biexp( Indo.1$time, Phi ), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(Indo.1), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSbiexp( Indo.1$time, A1, lrc1, A2, lrc2 ),
                        SS_biexp( Indo.1$time, Phi ), check.attributes = FALSE))
})

context("nlmer")
test_that("SSfol works", {
  Theoph.1 <- Theoph[ Theoph$Subject == 1, ]
  lKe <- -2.5; lKa <- 0.5; lCl <- -3
  Phi <- cbind(lKe, lKa, lCl)
  expect_true(all.equal(SSfol(Theoph.1$Dose, Theoph.1$Time, lKe, lKa, lCl),
                        SS_fol(Theoph.1$Dose, Theoph.1$Time, Phi), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(Theoph.1), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSfol(Theoph.1$Dose, Theoph.1$Time, lKe, lKa, lCl),
                        SS_fol(Theoph.1$Dose, Theoph.1$Time, Phi), check.attributes = FALSE))
})

context("nlmer")
test_that("SSfpl works", {
  Chick.1 <- ChickWeight[ChickWeight$Chick == 1, ]
  A <- 13; B <- 368; xmid <- 14; scal <- 6
  Phi <- cbind(A, B, xmid, log(scal))
  expect_true(all.equal(SSfpl(Chick.1$Time, A, B, xmid, scal),
                        SS_fpl(Chick.1$Time, Phi), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(Chick.1), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSfpl(Chick.1$Time, A, B, xmid, scal),
                        SS_fpl(Chick.1$Time, Phi), check.attributes = FALSE))
})

context("nlmer")
test_that("SSgompertz works", {
  DNase.1 <- subset(DNase, Run == 1)
  Asym <- 4.5; b2 <- 2.3; b3 <- 0.7
  Phi <- cbind(Asym, b2, b3)
  expect_true(all.equal(SSgompertz(log(DNase.1$conc), Asym, b2, b3),
                        SS_gompertz(log(DNase.1$conc), Phi), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(DNase.1), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSgompertz(log(DNase.1$conc), Asym, b2, b3),
                        SS_gompertz(log(DNase.1$conc), Phi), check.attributes = FALSE))
})

context("nlmer")
test_that("SSlogis works", {
  Chick.1 <- ChickWeight[ChickWeight$Chick == 1, ]
  Asym <- 368; xmid <- 14; scal <- 6
  Phi <- cbind(Asym, xmid, log(scal))
  expect_true(all.equal(SSlogis(Chick.1$Time, Asym, xmid, scal),
                        SS_logis(Chick.1$Time, Phi), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(Chick.1), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSlogis(Chick.1$Time, Asym, xmid, scal),
                        SS_logis(Chick.1$Time, Phi), check.attributes = FALSE))
})

context("nlmer")
test_that("SSmicmen works", {
  PurTrt <- Puromycin[ Puromycin$state == "treated", ]
  Vm <- 200; K <- 0.05
  Phi <- cbind(Vm, K)
  expect_true(all.equal(SSmicmen(PurTrt$conc, Vm, K),
                        SS_micmen(PurTrt$conc, Phi), check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(PurTrt), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSmicmen(PurTrt$conc, Vm, K),
                        SS_micmen(PurTrt$conc, Phi), check.attributes = FALSE))
})

context("nlmer")
test_that("SSweibull works", {
  Chick.6 <- subset(ChickWeight, (Chick == 6) & (Time > 0))
  Asym <- 160; Drop <- 115; lrc <- -5.5; pwr <- 2.5
  Phi <- cbind(Asym, Drop, lrc, pwr)
  expect_true(all.equal(SSweibull(Chick.6$Time, Asym, Drop, lrc, pwr) ,
                        SS_weibull(Chick.6$Time, Phi) , check.attributes = FALSE))
  Phi <- matrix(Phi, nrow = nrow(Chick.6), ncol = ncol(Phi), byrow = TRUE)
  expect_true(all.equal(SSweibull(Chick.6$Time, Asym, Drop, lrc, pwr) ,
                        SS_weibull(Chick.6$Time, Phi) , check.attributes = FALSE))
})

context("nlmer")
test_that("reshape works", {
  x <- as.double(1:10)
  expect_true(all(matrix(x, 5, 2) == reshape_vec(x, 5L, 2L)))
})

# betareg
links <- c("logit", "probit", "cloglog", "cauchit", "log")

context("betareg")
test_that("linkinv_beta returns expected results", {
  for (i in 1:length(links)) {
    eta <- -abs(rnorm(N))
    linkinv <- binomial(link = links[i])$linkinv
    expect_true(all.equal(linkinv(eta), 
                          linkinv_beta(eta, i)), info = links[i])
  }
})
context("betareg")
test_that("pw_beta and ll_beta_lp return expected results", {
  for (i in 1:length(links)) {
    eta <- -abs(rnorm(N))
    mu <- linkinv_beta(eta, i)
    dispersion <- 4/3
    linkinv <- binomial(link = links[i])$linkinv
    ll <- dbeta(1/3, mu*dispersion, (1-mu)*dispersion, log = TRUE)
    expect_true(all.equal(ll, pw_beta(rep(1/3,N) , eta, dispersion, i)), info = links[i])
  }
})

context("clogit")
test_that("ll_clogit_lp (which calls log_clogit_denom) returns the expected results", {
  data(infert)
  infert <- infert[order(infert$stratum, !infert$case),]
  betas <- c(spontaneous = 1.98587551667772, induced = 1.40901163187514)
  X <- model.matrix(case ~ spontaneous + induced - 1, data = infert)
  eta <- c(X %*% betas)
  y <- infert$case == 1
  s <- aggregate(y, by = list(infert$stratum), FUN = sum)$x
  obs <- aggregate(y, by = list(infert$stratum), FUN = length)$x
  ll <- ll_clogit_lp(eta0 = eta[!y], eta1 = eta[y], 
                     successes = s, failures = obs - s, observations = obs)
  expect_equal(-64.202236924431, ll)
})
