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
  if(!log) return(exp(log_f))
  else return(log_f)
}

loglog <- list(linkfun = qgumbel, linkinv = pgumbel, mu.eta = dgumbel, 
               valideta = function(eta) TRUE, name = "loglog")
class(loglog) <- "link-glm"

#' Fitting ordinal regression models via Stan
#'
#' Full Bayesian inference or optimization for ordinal regression models with
#' Gaussian, Student t, or Cauchy prior distributions for the coefficients.
#'
#' @export
#'
#' @param formula,data,weights,subset,na.action,contrasts,model,method 
#'   Same as in \code{\link[MASS]{polr}}.
#' @param Hess Same as in \code{\link[MASS]{polr}} but always taken to be
#'   \code{TRUE} and moreover ignored in the case of MCMC
#' @param ... Further arguments passed to \code{\link[rstan]{sampling}} (e.g.
#'   \code{iter}, \code{chains}, \code{refresh}, etc.) or 
#'   \code{\link[rstan]{optimizing}}
#' @param start Same as \code{\link[stats]{glm}}, but if not \code{NULL} also
#'   used as starting values for the MCMC. If \code{NULL} (the default), then
#'   \code{\link[rstan]{stan}} is initialized with \code{init = 'random'}.
#'   
#' @param prior Prior for parameters. Can be \code{NULL} to omit a prior
#'   and see \code{\link{priors}} otherwise.
#' @param prior_counts A numeric vector that must be positive but need not
#'   contain integers representing the prior count in each outcome. 
#'   Can be \code{NULL} to use a uniform Dirichlet prior.
#' @param prior_PD A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to draw from the prior predictive distribution instead of
#'   conditioning on the outcome
#' @param algorithm Character string (possibly abbreviated) among 
#'   \code{"sampling"} and \code{"optimizing"} indicating what estimation
#'   approach to use.
#'
#' @details The \code{stan_polr} function is similar in syntax to 
#'   \code{\link[MASS]{polr}} but rather than performing maximum likelihood 
#'   estimation of a proportional odds model, full Bayesian estimation is (by
#'   default) performed via Markov Chain Monte Carlo. The \code{stan_polr} 
#'   function calls the workhorse \code{stan_polr.fit} function, but it is 
#'   possible to call the latter directly.
#' 
#' @return A named list containing some components
#' @examples 
#' \dontrun{
#' # coming soon
#' }

stan_polr <- function (formula, data, weights, start, ..., subset, 
                       na.action = getOption("na.action", "na.omit"), 
                       contrasts = NULL, Hess = FALSE, model = TRUE, 
                       method = c("logistic", "probit", "loglog", "cloglog", "cauchit"),
                       prior = R2(stop("'location' must be specified")), 
                       prior_counts = NULL, prior_PD = FALSE, 
                       algorithm = c("sampling", "optimizing")) {
  
  # parse it like polr does in the MASS package
  m <- match.call(expand.dots = FALSE)
  method <- match.arg(method)
  if (is.matrix(eval.parent(m$data))) 
    m$data <- as.data.frame(data)
  m$start <- m$Hess <- m$method <- m$model <- m$... <- m$prior <- m$prior_counts <- 
    m$prior_PD <- m$algorithm <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
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
  else warning("an intercept is needed and assumed")
  K <- ncol(x)
  wt <- model.weights(m)
  if (!length(wt)) wt <- rep(1, n)
  offset <- model.offset(m)
  if (length(offset) <= 1L) offset <- rep(0, n)
  y <- model.response(m)
  if (!is.factor(y)) stop("response must be a factor")
  lev <- levels(y)
  llev <- length(lev)
  if (llev < 2L) stop("response must have 2 or more levels")
  # y <- unclass(y)
  q <- llev - 1L
  if (missing(start)) start <- NULL
  else if (length(start) != pc + q) 
    stop("'start' is not of the correct length")
  
  algorithm <- match.arg(algorithm)  
  stanfit <- stan_polr.fit(x, y, wt, start, offset, method, 
                           prior = prior, prior_counts = prior_counts,
                           prior_PD = prior_PD, algorithm = algorithm, ...)

  # list of all the arguments and their values including any defaults (match.call
  # doesn't include defaults)
  all_args <- mget(names(formals()), sys.frame(sys.nframe()))
  prior.info <- all_args[grep("prior", names(all_args), fixed = TRUE)]
  call <- match.call(expand.dots = TRUE)

  if (method == "logistic") linkinv <- plogis
  else if (method == "loglog") linkinv <- pgumbel
  else linkinv <- make.link(method)$linkinv

  rank <- qr(x, tol = .Machine$double.eps, LAPACK = TRUE)$rank
  df.residual <- n - sum(wt == 0) - rank
  
  if (llev == 2) { # actually a Bernoulli model
    if      (method == "logistic") family <- binomial(link = "logit")
    else if (method == "loglog")   family <- binomial(loglog)
    else                           family <- binomial(link = method)
    
    fit <- nlist(stanfit, family, formula, offset, weights = wt,
                 x = cbind("(Intercept)" = 1, x), y = as.integer(y == lev[2]), 
                 data, prior.info, call, terms = Terms, model = m,
                 na.action = attr(m, "na.action"), 
                 contrasts = attr(x, "contrasts"))
    out <- stanreg(fit)
    if (!model) out$model <- NULL
    class(out) <- c("stanreg", "polr")
    return(out)
  }
  else if (algorithm == "optimizing") {
    L <- t(chol(stanfit$cov.scaled))
    k <- nrow(L)
    unconstrained <- stanfit$par[1:k] + L %*% matrix(rnorm(4000 * k), k)
    stanmat <- t(apply(unconstrained, 2, FUN = function(u)
      unlist(constrain_pars(stanfit$stanfit, u))))
    stan_summary <- cbind(Estimate = stanfit$par, 
                          "Std. Error" = apply(stanmat, 2, sd),
                          t(apply(stanmat, 2, quantile,
                                  probs = c(0.025, .975))))
    covmat <- cov(stanmat)
    coefs <- stanfit$par[colnames(x)]
    eta <- linear_predictor(coefs, x, offset)
    mu <- linkinv(eta)
    residuals <- NA_real_
    names(eta) <- names(mu) <- rownames(x)
    
    llargs <- nlist(family = method, x, y, weights = wt, offset, 
                    beta = stanmat[,grep("^beta[[:digit:]]+$", colnames(stanmat)),drop = FALSE], 
                    zeta = stanmat[,grep("^zeta[[:digit:]]+$", colnames(stanmat)),drop=FALSE])
    log_lik <- do.call("pw_log_lik", llargs)
  }
  else {
    stanmat <- as.matrix(stanfit)
    coefs <- colMeans(stanmat[,1:K,drop=FALSE])
    eta <- linear_predictor(coefs, x, offset)
    mu <- linkinv(eta)
    
    means <- rstan::get_posterior_mean(stanfit)
    residuals <- means[grep("^residuals\\[[[:digit:]]+\\]$", rownames(means)),ncol(means)]
    names(residuals) <- names(eta) <- names(mu) <- rownames(x)
    
    llargs <- nlist(family = method, x, y, weights = wt, offset, 
                    beta = stanmat[,colnames(x),drop = FALSE], 
                    zeta = stanmat[,grep("|", colnames(stanmat), fixed = TRUE, value = TRUE),drop=FALSE])
    log_lik <- do.call("pw_log_lik", llargs)
    
    levs <- c(0.5, 0.8, 0.95, 0.99)
    qq <- (1 - levs) / 2
    probs <- sort(c(0.5, qq, 1 - qq))
    stan_summary <- rstan::summary(stanfit, probs = probs, digits = 10)$summary
  }
  
  out <- list(coefficients = coefs, fitted.values = mu, linear.predictors = eta,
              residuals = residuals, df.residual = df.residual, covmat = cov(stanmat),
              y = y, x = x, model = if(model) m, data = data, rank = rank,
              offset = offset, weights = wt, prior.weights = wt,
              method = method, contrasts = contrasts, na.action = na.action,
              call = call, formula = formula,
              terms = Terms, prior.info = prior.info, log_lik = log_lik,
              algorithm = algorithm,
              stan_summary = stan_summary, 
              stanfit = if (algorithm == "optimizing") stanfit$stanfit else stanfit)
  class(out) <- c("stanreg", "polr")
  return(out)
}
