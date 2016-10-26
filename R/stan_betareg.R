# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
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

#' Bayesian beta regression models via Stan
#'
#' Beta regression modeling with optional prior distributions for 
#' the coefficients, intercept, and dispersion parameter.
#'
#' @export
#' @templateVar armRef (Ch. 3-6)
#' @templateVar pkg stats
#' @templateVar pkgfun betareg
#' @templateVar sameargs model,offset,weights 
#' @templateVar rareargs na.action,contrasts
#' @templateVar fun stan_betareg
#' @templateVar fitfun stan_betareg.fit
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template see-also
#' @template args-formula-data-subset
#' @template args-same-as
#' @template args-same-as-rarely
#' @template args-x-y
#' @template args-dots
#' @template args-priors
#' 
#' @param prior_z See \code{prior}.
#' @param prior_intercept_z See \code{prior_intercept}. 
#' 
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template args-sparse
#' @template reference-gelman-hill
#' 
#' @param link Same as \code{\link[betareg]{betareg}}.
#' @param link.phi Same as \code{\link[betareg]{betareg}}
#' 
#' @details The \code{stan_betareg} function is similar in syntax to 
#'   \code{\link[betareg]{betareg}} but rather than performing maximum likelihood 
#'   estimation of generalized linear models, full Bayesian estimation is 
#'   performed (if \code{algorithm} is \code{"sampling"}) via MCMC. The Bayesian
#'   model adds independent priors on the coefficients of the beta regression model. The 
#'   \code{stan_betareg} function calls the workhorse \code{stan_betareg.fit} function, 
#'   but it is also possible to call the latter directly.
#'   
#' @seealso The various vignettes for \code{stan_betareg}.
#' 
#' @examples 
#' \donttest{
#' ### Simulated data
#' fam <- binomial(link = "logit")
#' dat <- list()
#' dat$N <- 200
#' dat$x <- rnorm(dat$N, 2, 1)
#' dat$z <- rnorm(dat$N, 2, 1)
#' dat$mu <- fam$linkinv(1 + 0.2*dat$x)
#' dat$phi <- exp(1.5 + 0.4*dat$z)
#' dat$y <- rbeta(dat$N, dat$mu * dat$phi, (1 - dat$mu) * dat$phi)
#' hist(dat$y, col = "dark grey", border = F, xlim = c(0,1))
#' fake_dat <- data.frame(dat$y, dat$x, dat$z)
#' colnames(fake_dat) <- c("y", "x", "z")
#' 
#' fit <- stan_betareg(y ~ x | z, data = fake_dat, link = "logit",
#'                     link.phi = "log", algorithm = "sampling")
#' print(fit, digits = 2)
#' pp_check(fit)
#' }

stan_betareg <- function (formula, data, subset, na.action, weights, offset,
                          link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                          link.phi = c("log", "identity", "sqrt"), model = TRUE, y = TRUE, x = FALSE, ...,
                          prior = normal(), prior_intercept = normal(),
                          prior_z = normal(), prior_intercept_z = normal(),
                          prior_ops = prior_options(), prior_PD = FALSE, 
                          algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                          adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  mc <- match.call(expand.dots = FALSE)
  mc$model <- mc$y <- mc$x <- TRUE
  
  # NULLify any Stan specific arguments in mc now
  mc$prior <- mc$prior_intercept <- mc$prior_ops <- mc$prior_PD <- mc$algorithm <-
    mc$adapt_delta <- mc$QR <- mc$sparse <- NULL
  
  mc$drop.unused.levels <- TRUE
  mc[[1L]] <- quote(betareg::betareg)
  
  if (!requireNamespace("betareg")) stop("the betareg package is needed by 'stan_betareg'")
  mc$control <- betareg::betareg.control(maxit = 0, fsmaxit = 0)
  br <- suppressWarnings(eval(mc, parent.frame()))
  mf <- check_constant_vars(br$model)
  mt <- br$terms
  Y <- array1D_check(model.response(mf, type = "any"))
  X <- model.matrix(br)
  Z <- model.matrix(br, model = "precision")
  #if(ncol(Z) == 1 && all(Z == 1)) Z <- NULL
  
  weights <- validate_weights(as.vector(model.weights(mf)))
  offset <- validate_offset(as.vector(model.offset(mf)), y = Y)
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  
  # pass existence of declaration of linear predictor of the dispertion parameter

  Z_true <- length(grep("\\|", all.names(formula)))
  # Z_true <- all.names(formula)[3]
  # if (Z_true=="|") {
  #   Z_true <- 1
  # }
  # else {
  #   Z_true <- 0
  # }

  
  # pass the prior information to stan_betareg.fit()
  stanfit <- stan_betareg.fit(x = X, y = Y, z = Z, weights = NULL, offset = NULL, 
                              link = link, link.phi = link.phi, ..., prior = prior, 
                              prior_intercept = prior_intercept, 
                              prior_z = prior_z, prior_intercept_z = prior_intercept_z,
                              prior_ops = prior_ops,
                              prior_PD = prior_PD, algorithm = algorithm, 
                              adapt_delta = adapt_delta, QR = QR, sparse = FALSE, Z_true = Z_true)
  algorithm <- match.arg(algorithm)
  link <- match.arg(link)
  link_phi <- match.arg(link.phi)
  fit <- nlist(stanfit, family = beta_fam(link), family_phi = beta_phi_fam(link_phi), formula, offset = NULL, 
               weights = NULL, x = X, y = Y, z = Z, # need Z if it has 2+ columns
               data, prior.info = get_prior_info(call, formals()), 
               call = match.call(), terms = mt, model = mf, 
               algorithm, na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"))
  out <- stanreg(fit)
  class(out) <- c("stanreg", "betareg")
  # out$xlevels <- .getXlevels(mt, mf)
  if (!x) 
    out$x <- NULL
  if (!y) 
    out$y <- NULL
  if (!model)
    out$model <- NULL
  
  return(out) 
}

beta_fam <- function(link = "logit") { # change function name to beta_mu_fam
  stopifnot(is.character(link))
  if (link == "loglog") {
    out <- binomial("cloglog")
    out$linkinv <- function(eta) 1 - pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
    out$linkfun <- function(mu) log(-log(mu))
  }
  else out <- binomial(link)
  out$family <- "beta"
  out$variance <- function(mu, phi) mu * (1 - mu) / (phi + 1)
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}

beta_phi_fam <- function(link = "log") {
  stopifnot(is.character(link))
  out <- poisson(link)
  out$family <- "beta_phi"
  out$variance <- function(mu, phi) mu * (1 - mu) / (phi + 1)
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
  return(out)  
}