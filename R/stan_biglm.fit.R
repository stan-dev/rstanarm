# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016 Trustees of Columbia University
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

#' @rdname stan_biglm
#' @export
#' @param b A numeric vector of OLS coefficients --- excluding the intercept --- 
#'   or a list of such vectors with one list element for each group
#' @param R A square upper-triangular matrix from the QR decomposition of the 
#'   centered design matrix or a list of such matrices with one list element for
#'   each group
#' @param SSR A numeric scalar indicating the sum-of-squared residuals for OLS
#'   or a numeric vector of sums-of-squared residuals with one element for each
#'   group
#' @param N A integer scalar indicating the number of included observations or
#'   an integer array of sample sizes with one element for each group
#' @template args-kappa_mean   
#' @examples
#' # create inputs
#' ols <- lm(mpg ~ wt + qsec + am - 1, # next line is critical for centering
#'           data = as.data.frame(scale(mtcars, scale = FALSE)))
#' b <- coef(ols)
#' R <- qr.R(ols$qr)
#' SSR <- crossprod(ols$residuals)[1]
#' N <- length(ols$fitted.values)
#' xbar <- colMeans(model.matrix(ols))
#' y <- mtcars$mpg
#' ybar <- mean(y)
#' s_y <- sd(y)
#' post <- stan_biglm.fit(b, R, SSR, N, xbar, ybar, s_y, prior = R2(.75),
#'                        # the next line is only to make the example go fast
#'                        chains = 1, iter = 1000, seed = 12345)
#' cbind(lm = b, stan_lm = rstan::get_posterior_mean(post)[14:16]) # shrunk
stan_biglm.fit <- function(b, R, SSR, N, xbar, ybar, s_y, has_intercept = TRUE, ...,
                           prior = R2(stop("'location' must be specified")), 
                           prior_intercept = NULL, kappa_mean = 1, prior_PD = FALSE,
                           algorithm = c("sampling", "meanfield", "fullrank"),
                           adapt_delta = NULL) {
  if (is.list(b)) {
    J <- length(b)
    K <- ncol(R[[1]])
    cn <- names(xbar[[1]])
    if (is.null(cn)) cn <- names(b[[1]])
    I <- diag(K)
    R_inv <- sapply(R, simplify = FALSE, FUN = backsolve, x = I)
    center_y <- 0.0
    xbarR_inv <- Rb <- vector("list", J)
    names(xbarR_inv) <- names(Rb) <- names(b)
    for (j in 1:J) {
      xbarR_inv[[j]] <- c(xbar[[j]] %*% R_inv[[j]])
      Rb[[j]] <- c(R[[j]] %*% b[[j]])
    }
  }
  else {
    J <- 1L
    N <- array(N, c(J))
    K <- ncol(R)
    cn <- names(xbar)
    if (is.null(cn)) cn <- names(b)
    I <- diag(K)
    R_inv <- backsolve(R, I)
    JK <- c(J, K)
    xbarR_inv <- array(c(xbar %*% R_inv), JK)
    Rb <- array(R %*% b, JK)
    SSR <- array(SSR, J)
    s_y <- array(s_y, J)
    center_y <- if (isTRUE(all.equal(matrix(0, J, K), xbar))) ybar else 0.0
    ybar <- array(ybar, J)
    dim(R_inv) <- c(J, dim(R_inv))
  }
  
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
  
  # initial values
  if (J == 1) R2 <- array(1 - SSR[1] / ((N - 1) * s_y^2), J)
  else R2 <- 1 - SSR / ((N - 1) * s_y^2)
  log_omega <- array(0, ifelse(prior_PD == 0, J, 0))
  init_fun <- function(chain_id) {
    out <- list(R2 = R2, log_omega = log_omega)
    if (has_intercept == 0L) out$z_alpha <- double()
    if (J == 1) {
      out$SMC <- double()
      out$mu <- array(0, c(0, K))
      out$kappa <- double()
    }
    return(out)
  }
  stanfit <- stanmodels$lm
  standata <- nlist(K, has_intercept, prior_dist,
                    prior_dist_for_intercept, 
                    prior_mean_for_intercept,
                    prior_scale_for_intercept, kappa_mean,
                    prior_PD, eta, J, N, xbarR_inv,
                    ybar, center_y, s_Y = s_y, Rb, SSR, R_inv)
  pars <- c(if (has_intercept) "alpha", "beta", "sigma", 
            if (prior_PD == 0) "log_omega", "R2", "mean_PPD")
  algorithm <- match.arg(algorithm)
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
  if (J == 1) new_names <- c(if (has_intercept) "(Intercept)", cn, "sigma", 
                             if (prior_PD == 0) "log-fit_ratio", 
                             "R2", "mean_PPD", "log-posterior")
  else {
    group_names <- names(b)
    new_names <- c(if (has_intercept) paste0("(Intercept):", group_names),
                   t(sapply(group_names, FUN = function(g) paste0(cn, ":", g))),
                   paste0("sigma:", group_names),
                   if (prior_PD == 0) paste0("log-fit_ratio:", group_names),
                   paste0("R2:", group_names), paste0("mean_PPD:", group_names),
                   "log-posterior")
  }
  stanfit@sim$fnames_oi <- new_names
  return(stanfit)
}
