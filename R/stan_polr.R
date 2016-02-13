# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright 1994-2013 William N. Venables and Brian D. Ripley
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

pgumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE) {
  q <- (q - loc)/scale
  p <- exp(-exp(-q))
  if (!lower.tail) 
    1 - p
  else 
    p
}

qgumbel <- function(p, loc = 0, scale = 1) {
  loc - scale * log(-log(p))
}

dgumbel <- function(x, loc = 0, scale = 1, log = FALSE) {
  z <- (x - loc) / scale
  log_f <- -(z + exp(-z))
  if (!log) 
    exp(log_f)
  else 
    log_f
}

loglog <- list(linkfun = qgumbel, linkinv = pgumbel, mu.eta = dgumbel, 
               valideta = function(eta) TRUE, name = "loglog")
class(loglog) <- "link-glm"

#' Bayesian ordinal regression models via Stan
#'
#' Bayesian inference for ordinal (or binary) regression models under
#' a proportional odds assumption.
#'
#' @export
#' @templateVar fun stan_polr
#' @templateVar fitfun stan_polr.fit
#' @templateVar pkg MASS
#' @templateVar pkgfun polr
#' @templateVar rareargs weights,na.action,contrasts,model
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template see-also
#' @template args-formula-data-subset
#' @template args-same-as-rarely
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-dots
#' @template args-adapt_delta
#'
#' @param method One of 'logistic', 'probit', 'loglog', 'cloglog' or 'cauchit',
#'   but can be abbreviated. See \code{\link[MASS]{polr}} for more details.
#' @param prior Prior for coefficients. Should be a call to \code{\link{R2}} 
#'   to specify the prior location of the \eqn{R^2} but can be \code{NULL}
#'   to indicate a standard uniform prior. See \code{\link{priors}}.
#' @param prior_counts A call to \code{\link{dirichlet}} to specify the 
#'   prior counts of the outcome when the predictors are at their sample
#'   means.
#' @param shape Either \code{NULL} or a positive scalar that is interpreted
#'   as the shape parameter for a \code{\link[stats]{GammaDist}}ribution on
#'   the exponent applied to the probability of success when there are only
#'   two outcome categories. If \code{NULL}, which is the default, then the
#'   exponent is taken to be fixed at \eqn{1}.
#' @param rate Either \code{NULL} or a positive scalar that is interpreted
#'   as the rate parameter for a \code{\link[stats]{GammaDist}}ribution on
#'   the exponent applied to the probability of success when there are only
#'   two outcome categories. If \code{NULL}, which is the default, then the
#'   exponent is taken to be fixed at \eqn{1}.   
#'
#' @details The \code{stan_polr} function is similar in syntax to 
#'   \code{\link[MASS]{polr}} but rather than performing maximum likelihood 
#'   estimation of a proportional odds model, Bayesian estimation is performed
#'   (if \code{algorithm = "sampling"}) via MCMC. The \code{stan_polr} 
#'   function calls the workhorse \code{stan_polr.fit} function, but it is 
#'   possible to call the latter directly.
#'   
#'   As for \code{\link{stan_lm}}, it is necessary to specify the prior 
#'   location of \eqn{R^2}. In this case, the \eqn{R^2} pertains to the
#'   proportion of variance in the latent variable (which is discretized
#'   by the cutpoints) attributable to the predictors in the model. 
#'   
#'   Prior beliefs about the cutpoints are governed by prior beliefs about the
#'   outcome when the predictors are at their sample means. Both of these
#'   are explained in the help page on \code{\link{priors}} and in the 
#'   \pkg{rstanarm} vignettes.
#'  
#'   Unlike \code{\link[MASS]{polr}}, \code{stan_polr} also allows the "ordinal"
#'   outcome to contain only two levels, in which case the likelihood is the
#'   same by default as for \code{\link{stan_glm}} with \code{family = binomial}
#'   but the prior on the coefficients is different. However, \code{stan_polr}
#'   allows the user to specify the \code{shape} and \code{rate} hyperparameters,
#'   in which case the probability of success is defined as the logistic CDF of
#'   the linear predictor, raised to the power of \code{alpha} where \code{alpha}
#'   has a gamma prior with the specified \code{shape} and \code{rate}. This
#'   likelihood is called \dQuote{scobit} by Nagler (1994) because if \code{alpha}
#'   is not equal to \eqn{1}, then the relationship between the linear predictor
#'   and the probability of success is skewed. If \code{shape} or \code{rate} is
#'   \code{NULL}, then \code{alpha} is assumed to be fixed to \eqn{1}. 
#'   
#'   Otherwise, it is usually advisible to set \code{shape} and \code{rate} to
#'   the same number so that the expected value of \code{alpha} is \eqn{1} while
#'   leaving open the possibility that \code{alpha} may depart from \eqn{1} a
#'   little bit. It is often necessary to have a lot of data in order to estimate
#'   \code{alpha} with much precision and always necessary to inspect the
#'   Pareto shape parameters calculated by \code{\link{loo}} to see if the 
#'   results are particularly sensitive to individual observations.
#'   
#'   Users should think carefully about how the outcome is coded when using 
#'   a scobit-type model. When \code{alpha} is not \eqn{1}, the asymmetry
#'   implies that the probability of success is most sensitive to the predictors
#'   when the probability of success is less than \eqn{0.63}. Reversing the
#'   coding of the successes and failures allows the predictors to have the
#'   greatest impact when the probability of failure is less than \eqn{0.63}.
#'   Also, the gamma prior on \code{alpha} is positively skewed, but you
#'   can reverse the coding of the successes and failures to circumvent this
#'   property.
#'   
#' @references 
#' Nagler, J., (1994). Scobit: An Alternative Estimator to Logit and Probit.
#' \emph{American Journal of Political Science}. 230 -- 255.
#' 
#' @seealso The vignette for \code{stan_polr}.
#'
#' @examples
#' if (!grepl("^sparc",  R.version$platform))
#' stan_polr(tobgp ~ agegp, data = esoph, method = "probit",
#'           prior = R2(0.2, "mean"), init_r = 0.1, seed = 12345,
#'           algorithm = "fullrank") # for speed only
#' 
stan_polr <- function(formula, data, weights, ..., subset, 
                      na.action = getOption("na.action", "na.omit"), 
                      contrasts = NULL, model = TRUE, 
                      method = c("logistic", "probit", "loglog", "cloglog", 
                                 "cauchit"),
                      prior = R2(stop("'location' must be specified")), 
                      prior_counts = dirichlet(1), shape = NULL, rate = NULL, 
                      prior_PD = FALSE, 
                      algorithm = c("sampling", "meanfield", "fullrank"),
                      adapt_delta = NULL) {
  
  algorithm <- match.arg(algorithm)
  call <- match.call(expand.dots = TRUE)
  m <- match.call(expand.dots = FALSE)
  method <- match.arg(method)
  if (is.matrix(eval.parent(m$data))) 
    m$data <- as.data.frame(data)
  m$method <- m$model <- m$... <- m$prior <- m$prior_counts <- 
    m$prior_PD <- m$algorithm <- m$adapt_delta <- m$shape <- m$rate <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  m <- check_constant_vars(m)
  Terms <- attr(m, "terms")
  x <- model.matrix(Terms, m, contrasts)
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)
  n <- nrow(x)
  pc <- ncol(x)
  cons <- attr(x, "contrasts")
  if (xint > 0L) {
    x <- x[, -xint, drop = FALSE]
    pc <- pc - 1L
  } else {
    warning("An intercept is needed and assumed.")
  }
  K <- ncol(x)
  wt <- model.weights(m)
  if (!length(wt)) 
    wt <- rep(1, n)
  offset <- model.offset(m)
  if (length(offset) <= 1L) 
    offset <- rep(0, n)
  y <- model.response(m)
  if (!is.factor(y)) 
    stop("Response variable must be a factor.")
  lev <- levels(y)
  llev <- length(lev)
  if (llev < 2L) 
    stop("Response variable must have 2 or more levels.")
  # y <- unclass(y)
  q <- llev - 1L

  stanfit <- stan_polr.fit(x, y, wt, offset, method, 
                           prior = prior, prior_counts = prior_counts,
                           shape = shape, rate = rate,
                           prior_PD = prior_PD, algorithm = algorithm, 
                           adapt_delta = adapt_delta, ...)
  
  prior.info <- get_prior_info(call, formals())
  inverse_link <- linkinv(method)
  
  if (llev == 2L) { # actually a Bernoulli model
    if (method == "logistic") {
      family <- binomial(link = "logit")
    } else if (method == "loglog") {
      family <- binomial(loglog)
    } else {
      family <- binomial(link = method)
    }
    
    fit <- nlist(stanfit, family, formula, offset, weights = wt,
                 x = cbind("(Intercept)" = 1, x), y = as.integer(y == lev[2]), 
                 data, prior.info, call, terms = Terms, model = m,
                 algorithm, na.action = attr(m, "na.action"), 
                 contrasts = attr(x, "contrasts"))
    out <- stanreg(fit)
    if (!model) 
      out$model <- NULL
    class(out) <- c("stanreg", "polr")
    return(out)
  } else {
    # more than 2 outcome levels
    K2 <- K + llev - 1 # number of coefficients + number of cutpoints
    stanmat <- as.matrix(stanfit)[, 1:K2, drop = FALSE] 
    covmat <- cov(stanmat)
    coefs <- apply(stanmat[, 1:K, drop = FALSE], 2, median)
    ses <- apply(stanmat[, 1:K, drop = FALSE], 2, mad)
    zeta <- apply(stanmat[, (K+1):K2, drop = FALSE], 2, mad)
    eta <- linear_predictor(coefs, x, offset)
    mu <- inverse_link(eta)
    
    means <- rstan::get_posterior_mean(stanfit)
    residuals <- means[grep("^residuals", rownames(means)), ncol(means)]
    names(residuals) <- names(eta) <- names(mu) <- rownames(x)

    levs <- c(0.5, 0.8, 0.95, 0.99)
    qq <- (1 - levs) / 2
    probs <- sort(c(0.5, qq, 1 - qq))
    stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
    if (algorithm == "sampling") 
      check_rhats(stan_summary[, "Rhat"])
  }
  
  out <- nlist(coefficients = coefs, ses, zeta, residuals,
              fitted.values = mu, linear.predictors = eta, covmat,
              y, x, model = if (model) m, 
              family = method, data, offset, weights = wt, prior.weights = wt,
              method, contrasts, na.action,
              call, formula, terms = Terms, prior.info,
              algorithm, stan_summary, stanfit)
  class(out) <- c("stanreg", "polr")
  
  return(out)
}
