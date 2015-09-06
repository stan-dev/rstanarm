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

#' @rdname stan_lm
#' @export
stan_lm.wfit <- function(x, y, w, offset = NULL, singular.ok = TRUE, ...,
                         prior = R2(stop("'location' must be specified")), 
                         prior_PD = FALSE, 
                         algorithm = c("sampling", "optimizing", "meanfield",
                                       "fullrank")) {
  if (NCOL(y) > 1) stop("multivariate responses not supported yet")
  if (colnames(x)[1] == "(Intercept)") {
    has_intercept <- 1L
    x <- x[,-1,drop=FALSE]
  }
  else has_intercept <- 0L
  ols <- lsfit(x, y, w, has_intercept == 1L)
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
  if (isTRUE(all.equal(matrix(0, J, K), xbar))) center_y <- mean(y)
  else center_y <- 0
  ybar <- array(mean(y), J)

  eta <- prior$eta <- make_eta(prior$location, prior$what, K = K)
  
  # initial values
  L <- t(chol(cor(x)))
  R2 <- array(1 - SSR[1] / ((N - 1) * var(y)), J)
  log_omega <- array(0, ifelse(prior_PD == 0, J, 0))
  init_fun <- function(chain_id) {
    return(list(L = L, R2 = R2, log_omega = log_omega))
  }
  algorithm <- match.arg(algorithm)
  stanfit <- get("stanfit_lm")
  standata <- nlist(K, has_intercept, prior_PD, J, N, xbar, s_X, XtX, 
                    ybar, center_y, s_Y, b, SSR, eta)
  pars <- c(if (has_intercept) "alpha", "beta", "sigma", 
            if (prior_PD == 0) "log_omega", "mean_PPD")
  if (algorithm == "optimizing") {
    stop("'optimizing' not supported for this model")
  }
  else if (algorithm %in% c("meanfield", "fullrank")) {
    stanfit <- rstan::vb(stanfit, data = standata, pars = pars, ...)
  }
  else if ("control" %in% names(list(...))) {
    stanfit <- rstan::sampling(stanfit, data = standata, pars = pars, 
                               init = init_fun, show_messages = FALSE, ...)
  }
  else stanfit <- rstan::sampling(stanfit, data = standata, pars = pars, init = init_fun,
                                  control = list(adapt_delta = 0.95, max_treedepth = 15), 
                                  show_messages = FALSE, ...)
  parameters <- dimnames(stanfit)$parameters
  new_names <- c(if (has_intercept) "(Intercept)", colnames(x), "sigma", 
                 if (prior_PD == 0) "log-fit_ratio", "mean_PPD", "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  return(stanfit)
}

#' @rdname stan_lm
#' @export
stan_lm.fit <- function(x, y, offset = NULL, singular.ok = TRUE, ...,
                        prior = R2(stop("'location' must be specified")), 
                        prior_PD = FALSE, 
                        algorithm = c("sampling", "optimizing", "meanfield",
                                      "fullrank")) {

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("stan_lm.wfit")
  mf$w <- as.name("NULL")
  return(eval(mf, parent.frame()))
}
