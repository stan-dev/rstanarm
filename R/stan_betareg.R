# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Beta regression modeling with optional prior distributions for the 
#' coefficients, intercept, and auxiliary parameter \code{phi} (if applicable).
#'
#' @export
#' @templateVar armRef (Ch. 3-6)
#' @templateVar pkg betareg
#' @templateVar pkgfun betareg
#' @templateVar sameargs model,offset,weights 
#' @templateVar rareargs na.action
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
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' 
#' @param link Character specification of the link function used in the model 
#'   for mu (specified through \code{x}). Currently, "logit", "probit",
#'   "cloglog", "cauchit", "log", and "loglog" are supported.
#' @param link.phi If applicable, character specification of the link function 
#'   used in the model for \code{phi} (specified through \code{z}). Currently, 
#'   "identity", "log" (default), and "sqrt" are supported. Since the "sqrt"
#'   link function is known to be unstable, it is advisable to specify a
#'   different link function (or to model \code{phi} as a scalar parameter
#'   instead of via a linear predictor by excluding \code{z} from the
#'   \code{formula} and excluding \code{link.phi}).
#' @param prior_z Prior distribution for the coefficients in the model for 
#'   \code{phi} (if applicable). Same options as for \code{prior}.
#' @param prior_intercept_z Prior distribution for the intercept in the model 
#'   for \code{phi} (if applicable). Same options as for \code{prior_intercept}.
#' @param prior_phi The prior distribution for \code{phi} if it is \emph{not} 
#'   modeled as a function of predictors. If \code{z} variables are specified 
#'   then \code{prior_phi} is ignored and \code{prior_intercept_z} and 
#'   \code{prior_z} are used to specify the priors on the intercept and
#'   coefficients in the model for \code{phi}. When applicable, \code{prior_phi}
#'   can be a call to \code{exponential} to use an exponential distribution, or
#'   one of \code{normal}, \code{student_t} or \code{cauchy} to use half-normal,
#'   half-t, or half-Cauchy prior. See \code{\link{priors}} for details on these
#'   functions. To omit a prior ---i.e., to use a flat (improper) uniform
#'   prior--- set \code{prior_phi} to \code{NULL}.
#' 
#' @details The \code{stan_betareg} function is similar in syntax to 
#'   \code{\link[betareg]{betareg}} but rather than performing maximum 
#'   likelihood estimation, full Bayesian estimation is performed (if 
#'   \code{algorithm} is \code{"sampling"}) via MCMC. The Bayesian model adds 
#'   priors (independent by default) on the coefficients of the beta regression
#'   model. The \code{stan_betareg} function calls the workhorse
#'   \code{stan_betareg.fit} function, but it is also possible to call the
#'   latter directly.
#'   
#' @seealso The vignette for \code{stan_betareg}.
#'   \url{http://mc-stan.org/rstanarm/articles/}
#' 
#' @references Ferrari, SLP and Cribari-Neto, F (2004). Beta regression for 
#'   modeling rates and proportions. \emph{Journal of Applied Statistics}.
#'   31(7), 799--815.
#' 
#' @examples 
#' ### Simulated data
#' N <- 200
#' x <- rnorm(N, 2, 1)
#' z <- rnorm(N, 2, 1)
#' mu <- binomial(link = "logit")$linkinv(1 + 0.2*x)
#' phi <- exp(1.5 + 0.4*z)
#' y <- rbeta(N, mu * phi, (1 - mu) * phi)
#' hist(y, col = "dark grey", border = FALSE, xlim = c(0,1))
#' fake_dat <- data.frame(y, x, z)
#' 
#' fit <- stan_betareg(
#'   y ~ x | z, data = fake_dat, 
#'   link = "logit", 
#'   link.phi = "log", 
#'   algorithm = "optimizing" # just for speed of example
#'  ) 
#' print(fit, digits = 2)
#'
stan_betareg <-
  function(formula,
           data,
           subset,
           na.action,
           weights,
           offset,
           link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
           link.phi = NULL,
           model = TRUE,
           y = TRUE,
           x = FALSE,
           ...,
           prior = normal(autoscale=TRUE),
           prior_intercept = normal(autoscale=TRUE),
           prior_z = normal(autoscale=TRUE),
           prior_intercept_z = normal(autoscale=TRUE),
           prior_phi = exponential(autoscale=TRUE),
           prior_PD = FALSE,
           algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
           adapt_delta = NULL,
           QR = FALSE) {
    
    if (!requireNamespace("betareg", quietly = TRUE)) {
      stop("Please install the betareg package before using 'stan_betareg'.")
    }
    if (!has_outcome_variable(formula)) {
      stop("LHS of formula must be specified.")
    }
    
    mc <- match.call(expand.dots = FALSE)
    data <- validate_data(data, if_missing = environment(formula))
    mc$data <- data
    mc$model <- mc$y <- mc$x <- TRUE
    
    # NULLify any Stan specific arguments in mc
    mc$prior <- mc$prior_intercept <- mc$prior_PD <- mc$algorithm <-
      mc$adapt_delta <- mc$QR <- mc$sparse <- mc$prior_dispersion <- NULL
    
    mc$drop.unused.levels <- TRUE
    mc[[1L]] <- quote(betareg::betareg)
    mc$control <- betareg::betareg.control(maxit = 0, fsmaxit = 0)
    br <- suppressWarnings(eval(mc, parent.frame()))
    mf <- check_constant_vars(br$model)
    mt <- br$terms
    Y <- array1D_check(model.response(mf, type = "any"))
    X <- model.matrix(br)
    Z <- model.matrix(br, model = "precision")
    weights <- validate_weights(as.vector(model.weights(mf)))
    offset <- validate_offset(as.vector(model.offset(mf)), y = Y)
    
    # check if user specified matrix for precision model
    if (length(grep("\\|", all.names(formula))) == 0 && 
        is.null(link.phi))
      Z <- NULL
    
    algorithm <- match.arg(algorithm)
    link <- match.arg(link)
    link_phi <- match.arg(link.phi, c(NULL, "log", "identity", "sqrt"))
    
    stanfit <- 
      stan_betareg.fit(x = X, y = Y, z = Z, 
                       weights = weights, offset = offset,
                       link = link, link.phi = link.phi,
                       ...,
                       prior = prior, prior_z = prior_z,
                       prior_intercept = prior_intercept, 
                       prior_intercept_z = prior_intercept_z,
                       prior_phi = prior_phi, prior_PD = prior_PD,
                       algorithm = algorithm, adapt_delta = adapt_delta, 
                       QR = QR)
    if (algorithm != "optimizing" && !is(stanfit, "stanfit")) return(stanfit)
    if (is.null(link.phi) && is.null(Z))
      link_phi <- "identity"
    sel <- apply(X, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
    X <- X[ , !sel, drop = FALSE]
    if (!is.null(Z)) {
      sel <- apply(Z, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
      Z <- Z[ , !sel, drop = FALSE]
    }
    fit <- 
      nlist(stanfit, algorithm, data, offset, weights,
            x = X, y = Y, z = Z %ORifNULL% model.matrix(y ~ 1),
            family = beta_fam(link), family_phi = beta_phi_fam(link_phi),
            formula, model = mf, terms = mt, call = match.call(),
            na.action = attr(mf, "na.action"), contrasts = attr(X, "contrasts"), 
            stan_function = "stan_betareg")
    out <- stanreg(fit)
    if (algorithm == "optimizing") {
      out$log_p <- stanfit$log_p
      out$log_g <- stanfit$log_g
    }
    out$xlevels <- lapply(mf[,-1], FUN = function(x) {
      xlev <- if (is.factor(x) || is.character(x)) levels(x) else NULL
      xlev[!vapply(xlev, is.null, NA)]
    })
    out$levels <- br$levels
    if (!x)
      out$x <- NULL
    if (!y)
      out$y <- NULL
    if (!model)
      out$model <- NULL
    
    structure(out, class = c("stanreg", "betareg"))
  }


# internal ----------------------------------------------------------------
beta_fam <- function(link = "logit") {
  stopifnot(is.character(link))
  if (link == "loglog") {
    out <- binomial("cloglog")
    out$linkinv <- function(eta) {
      1 - pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), 
               .Machine$double.eps)
    }
    out$linkfun <- function(mu) log(-log(mu))
  } else {
    out <- binomial(link)
  }
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
}
