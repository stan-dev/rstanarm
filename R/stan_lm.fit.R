# This file is part of rstanarm.
# Copyright 2013 Stan Development Team
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

make_eta <- function(prior.R2, prior.what = c("mode", "mean", "median", "log"), K) {
  stopifnot(length(prior.R2) == 1, is.numeric(prior.R2))
  stopifnot(is.numeric(K), K > 0, K == as.integer(K))
  prior.what <- match.arg(prior.what)
  half_K <- K / 2
  if (prior.what == "mode") {
    stopifnot(prior.R2 > 0, prior.R2 <= 1)
    if (K <= 2) stop("mode of beta distribution does not exist when K <= 2 ",
                     "specify 'prior.what' as 'mean' or 'log' instead")
    eta <- (half_K - 1  - prior.R2 * half_K + prior.R2 * 2) / prior.R2
  }
  else if (prior.what == "mean") {
    stopifnot(prior.R2 > 0, prior.R2 <= 1)
    eta <- (half_K - prior.R2 * half_K) / prior.R2
  }
  else if (prior.what == "median") {
    stopifnot(prior.R2 > 0, prior.R2 <= 1)
    FUN <- function(eta) qbeta(0.5, half_K, qexp(eta)) - prior.R2
    eta <- qexp(uniroot(FUN, interval = 0:1)$root)
  }
  else { # prior.what == "log"
    stopifnot(prior.R2 < 0)
    FUN <- function(eta) digamma(half_K) - digamma(half_K + qexp(eta)) - prior.R2
    eta <- qexp(uniroot(FUN, interval = 0:1, 
                        f.lower = -prior.R2, f.upper = -.Machine$double.xmax)$root)
  }
  return(eta)  
}

#' @rdname stan_lm
#' @export
stan_lm.wfit <- function(x, y, w, offset = NULL, method = "qr", tol = 1e-07,
                    singular.ok = TRUE,
                    eta = NULL, prior.R2 = NULL, 
                    prior.what = c("mode", "mean", "median", "log"), ...) {
  
  if (colnames(x)[1] == "(Intercept)") {
    has_intercept <- 1L
    x <- x[,-1,drop=FALSE]
  }
  else has_intercept <- 0L
  ols <- lsfit(x, y, w, has_intercept == 1L, tol)
  if (!is.null(w)) x <- sqrt(w) * x

  J <- 1L
  N <- array(nrow(x), c(J))
  K <- ncol(x)
  if (K == 0) stop("'stan_lm.fit' is not suitable for estimating a mean ",
                   "use 'stan_glm.fit' with 'family = gaussian()' instead")
  b <- coef(ols)
  b[is.na(b)] <- 0.0
  if (has_intercept == 1L) b <- array(b[-1], c(J,K))
  else b <- array(b, c(J,K))
  
  SSR <- array(crossprod(residuals(ols))[1], J)
  
  s_X <- array(apply(x, 2, sd), c(J,K))
  xbar <- array(colMeans(x), c(J,K))
  x <- sweep(x, 2, xbar, FUN = "-")
  XtX <- crossprod(x)
  dim(XtX) <- c(J, K, K)
  s_Y <- array(sd(y), J)
  ybar <- array(mean(y), J)
  
  if (is.null(eta)) {
    if (is.null(prior.R2)) 
      stop("the 'prior.R2' argument must be specified if 'eta' is unspecified")
    eta <- make_eta(prior.R2, prior.what, K)
  }
  else if (!is.numeric(eta) || length(eta) != 1L || eta <= 0) {
    stop("'eta' must be a positive scalar")
  }
  
  stanfit <- get("stanfit_lm")
  standata <- nlist(K, has_intercept, J, N, xbar, s_X, XtX, ybar, s_Y, b, SSR, eta)
  pars <- c(if (has_intercept) "alpha", "beta", "sigma", "log_omega", "mean_PPD")
  stanfit <- rstan::sampling(stanfit, data = standata, pars = pars, ...)
  parameters <- dimnames(stanfit)$parameters
  new_names <- c(if (has_intercept) "(Intercept)", colnames(x), 
                 "sigma", "log-fit_ratio", "mean_PPD", "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  return(stanfit)
}

#' @rdname stan_lm
#' @export
stan_lm.fit <- function(x, y, offset = NULL, method = "qr", tol = 1e-07,
                         singular.ok = TRUE,
                         eta = NULL, prior.R2 = NULL, 
                         prior.what = c("mode", "mean", "median", "log"), ...) {

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("stan_lm.wfit")
  mf$w <- as.name("NULL")
  return(eval(mf, parent.frame()))
}
