# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015 Trustees of Columbia University
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

#' @rdname stan_lm
#' @export
stan_lm.wfit <- function(x, y, w, offset = NULL, singular.ok = TRUE, ...,
                         prior = R2(stop("'location' must be specified")), 
                         prior_intercept = NULL, prior_PD = FALSE, 
                         algorithm = c("sampling", "meanfield", "fullrank"),
                         adapt_delta = NULL) {
  algorithm <- match.arg(algorithm)
  if (NCOL(y) > 1) 
    stop("Multivariate responses not supported yet.")
  if (colnames(x)[1L] == "(Intercept)") {
    has_intercept <- 1L
    x <- x[, -1L, drop = FALSE]
    if (NCOL(x) == 0L)
      stop("'stan_lm' is not suitable for estimating a mean.",
           "\nUse 'stan_glm' with 'family = gaussian()' instead.", 
           call. = FALSE)
  } else {
    has_intercept <- 0L
  }
  if (nrow(x) < ncol(x))
    stop("stan_lm with more data points than predictors is not yet enabled.", 
         call. = FALSE)
  
  xbar <- colMeans(x)
  x <- sweep(x, 2L, xbar, FUN = "-")
  ybar <- mean(y)
  y <- y - ybar
  if(length(w) == 0) ols <- lm.fit(x, y)
  else ols <- lm.wfit(x, y, w)
  b <- coef(ols)
  NAs <- is.na(b)
  if (any(NAs) && !singular.ok) {
    x <- x[,!NAs, drop = FALSE]
    xbar <- xbar[!NAs]
    ols <- lsfit(x, y, w, intercept = FALSE)
    b <- coef(ols)
  }
  else b[NAs] <- 0.0
  
  if (!is.null(w)) 
    x <- sqrt(w) * x

  J <- 1L
  N <- array(nrow(x), c(J))
  K <- ncol(x)
  cn <- colnames(x)
  decomposition <- ols$qr
  Q <- qr.Q(decomposition)
  R <- qr.R(decomposition)
  R_inv <- qr.solve(decomposition, Q)
  x <- Q
  colnames(x) <- cn
  JK <- c(J, K)
  xbarR_inv <- array(c(xbar %*% R_inv), JK)
  Rb <- array(R %*% b, JK)

  SSR <- array(crossprod(residuals(ols))[1], J)
  s_Y <- array(sd(y), J)
  center_y <- if (isTRUE(all.equal(matrix(0, J, K), xbar))) ybar else 0
  ybar <- array(ybar, J)

  if (!length(prior)) {
    prior_dist <- 0L
    eta <- 0
  } else {
    prior_dist <- 1L
    eta <- prior$eta <- make_eta(prior$location, prior$what, K = K)
  }
  if (!length(prior_intercept)) {
    prior_dist_for_intercept <- 0L
    prior_mean_for_intercept <- 0
    prior_scale_for_intercept <- 0
  } else {
    if (!identical(prior_intercept$dist, "normal"))
      stop("'prior_intercept' must be 'NULL' or a call to 'normal'.")
    prior_dist_for_intercept <- 1L
    prior_mean_for_intercept <- prior_intercept$location
    prior_scale_for_intercept <- prior_intercept$scale
    if (is.null(prior_scale_for_intercept))
      prior_scale_for_intercept <- 0
  }
  dim(R_inv) <- c(J, dim(R_inv))
  
  # initial values
  R2 <- array(1 - SSR[1] / ((N - 1) * s_Y^2), J)
  log_omega <- array(0, ifelse(prior_PD == 0, J, 0))
  init_fun <- function(chain_id) {
    out <- list(R2 = R2, log_omega = log_omega)
    if (has_intercept == 0L) out$z_alpha <- double()
    return(out)
  }
  stanfit <- stanmodels$lm
  standata <- nlist(K, has_intercept, prior_dist,
                    prior_dist_for_intercept, 
                    prior_mean_for_intercept,
                    prior_scale_for_intercept,
                    prior_PD, eta, J, N, xbarR_inv,
                    ybar, center_y, s_Y, Rb, SSR, R_inv)
  pars <- c(if (has_intercept) "alpha", 
            "beta", 
            "sigma", 
            if (prior_PD == 0) "log_omega", 
            "R2", 
            "mean_PPD")
  if (algorithm %in% c("meanfield", "fullrank")) {
    stanfit <- rstan::vb(stanfit, data = standata, pars = pars,
                         algorithm = algorithm, ...)
  } else {
    sampling_args <- set_sampling_args(
      object = stanfit, 
      prior = prior,
      user_dots = list(...), 
      user_adapt_delta = adapt_delta, 
      init = init_fun, data = standata, pars = pars, show_messages = FALSE)
    stanfit <- do.call(sampling, sampling_args)
  }
  new_names <- c(if (has_intercept) "(Intercept)", 
                 colnames(x), 
                 "sigma", 
                 if (prior_PD == 0) "log-fit_ratio", 
                 "R2", 
                 "mean_PPD", 
                 "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  return(stanfit)
}

#' @rdname stan_lm
#' @export
stan_lm.fit <- function(x, y, offset = NULL, singular.ok = TRUE, ...,
                        prior = R2(stop("'location' must be specified")), 
                        prior_intercept = NULL, prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta = NULL) { # nocov start
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("stan_lm.wfit")
  mf$w <- as.name("NULL")
  eval(mf, parent.frame())
} # nocov end
