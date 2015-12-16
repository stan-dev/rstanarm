# This file is part of rstanarm.
# Copyright 1994-2013 William N. Venables and Brian D. Ripley
# Copyright 2015 Stan Development Team
# rstanarm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rstanarm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.

pgumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE) {
  q <- (q - loc)/scale
  p <- exp(-exp(-q))
  if (!lower.tail) 
    1 - p
  else p
}

qgumbel <- function(p, loc = 0, scale = 1) {
  loc - scale * log(-log(p))
}

dgumbel <- function(x, loc = 0, scale = 1, log = FALSE) {
  z <- (x - loc) / scale
  log_f <- -(z + exp(-z))
  if (!log) return(exp(log_f))
  else return(log_f)
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
#' @param prior Prior for coefficients. Can be \code{NULL} to omit a prior
#'   but otherwise must be a call to \code{\link{R2}} to specify the 
#'   prior location of the \eqn{R^2}. See \code{\link{priors}}.
#' @param prior_counts A call to \code{\link{dirichlet}} to specify the 
#'   prior counts of the outcome when the predictors are at their sample
#'   means.
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
#'   by the cutpoints) attributable to the predictors in the model. Prior
#'   beliefs about the cutpoints are governed by prior beliefs about the
#'   outcome when the predictors are at their sample means. Both of these
#'   are explained in the help page on \code{\link{priors}} and in the 
#'   \pkg{rstanarm} vignettes.
#'   
#'   Unlike \code{\link[MASS]{polr}}, \code{stan_polr} also allows the "ordinal"
#'   outcome to contain only two levels, in which case the likelihood is the
#'   same as for \code{\link{stan_glm}} with \code{family = binomial} but the
#'   prior on the coefficients is different.
#' 
#' @examples 
#' \dontrun{
#' stan_polr(tobgp ~ agegp, data = esoph, cores = 1,
#'           prior = R2(0.2, "mean"), init_r = 0.1, seed = 12345)
#' }
#' 
stan_polr <- function(formula, data, weights, ..., subset, 
                      na.action = getOption("na.action", "na.omit"), 
                      contrasts = NULL, model = TRUE, 
                      method = c("logistic", "probit", "loglog", "cloglog", 
                                 "cauchit"),
                      prior = R2(stop("'location' must be specified")), 
                      prior_counts = dirichlet(1), prior_PD = FALSE, 
                      algorithm = c("sampling", "meanfield", "fullrank"),
                      adapt_delta = NULL) {
  
  # parse it like polr does in the MASS package
  m <- match.call(expand.dots = FALSE)
  method <- match.arg(method)
  if (is.matrix(eval.parent(m$data))) 
    m$data <- as.data.frame(data)
  m$method <- m$model <- m$... <- m$prior <- m$prior_counts <- 
    m$prior_PD <- m$algorithm <- m$adapt_delta <- NULL
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
  }
  else warning("An intercept is needed and assumed.")
  K <- ncol(x)
  wt <- model.weights(m)
  if (!length(wt)) wt <- rep(1, n)
  offset <- model.offset(m)
  if (length(offset) <= 1L) offset <- rep(0, n)
  y <- model.response(m)
  if (!is.factor(y)) stop("Response variable must be a factor.")
  lev <- levels(y)
  llev <- length(lev)
  if (llev < 2L) stop("Response variable must have 2 or more levels.")
  # y <- unclass(y)
  q <- llev - 1L

  algorithm <- match.arg(algorithm)  
  stanfit <- stan_polr.fit(x, y, wt, offset, method, 
                           prior = prior, prior_counts = prior_counts,
                           prior_PD = prior_PD, algorithm = algorithm, 
                           adapt_delta = adapt_delta, ...)
  
  call <- match.call(expand.dots = TRUE)
  prior.info <- get_prior_info(call, formals())

  inverse_link <- linkinv(method)
  rank <- qr(x, tol = .Machine$double.eps, LAPACK = TRUE)$rank
  df.residual <- n - sum(wt == 0) - rank
  
  if (llev == 2) { # actually a Bernoulli model
    if (method == "logistic") family <- binomial(link = "logit")
    else if (method == "loglog") family <- binomial(loglog)
    else family <- binomial(link = method)
    
    fit <- nlist(stanfit, family, formula, offset, weights = wt,
                 x = cbind("(Intercept)" = 1, x), y = as.integer(y == lev[2]), 
                 data, prior.info, call, terms = Terms, model = m,
                 algorithm, na.action = attr(m, "na.action"), 
                 contrasts = attr(x, "contrasts"))
    out <- stanreg(fit)
    if (!model) out$model <- NULL
    class(out) <- c("stanreg", "polr")
    return(out)
  }
  else {
    K2 <- K + llev - 1 # number of coefficients + number of cutpoints
    stanmat <- as.matrix(stanfit)[, 1:K2, drop = FALSE] 
    covmat <- cov(stanmat)
    coefs <- apply(stanmat[, 1:K, drop = FALSE], 2, median)
    ses <- apply(stanmat[, 1:K, drop = FALSE], 2, mad)
    zeta <- apply(stanmat[, K:K2, drop = FALSE], 2, mad)
    eta <- linear_predictor(coefs, x, offset)
    mu <- inverse_link(eta)
    
    means <- rstan::get_posterior_mean(stanfit)
    residuals <- means[grep("^residuals", rownames(means)),ncol(means)]
    names(residuals) <- names(eta) <- names(mu) <- rownames(x)

    levs <- c(0.5, 0.8, 0.95, 0.99)
    qq <- (1 - levs) / 2
    probs <- sort(c(0.5, qq, 1 - qq))
    stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
  }
  
  out <- list(coefficients = coefs, ses = ses, zeta = zeta,
              fitted.values = mu, linear.predictors = eta,
              residuals = residuals, df.residual = df.residual, covmat = covmat,
              y = y, x = x, model = if (model) m, data = data, rank = rank,
              offset = offset, weights = wt, prior.weights = wt,
              method = method, contrasts = contrasts, na.action = na.action,
              call = call, formula = formula,
              terms = Terms, prior.info = prior.info,
              algorithm = algorithm,
              stan_summary = stan_summary, family = method,
              stanfit = stanfit)
  class(out) <- c("stanreg", "polr")
  return(out)
}

