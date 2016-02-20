# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

#' @rdname stan_polr
#' @export
#' @param x A design matrix.
#' @param y A response variable, which must be a (preferably ordered) factor.
#' @param wt A numeric vector (possibly \code{NULL}) of observation weights.
#' @param offset A numeric vector (possibly \code{NULL}) of offsets.
#' 
#' @importFrom utils head tail
stan_polr.fit <- function(x, y, wt = NULL, offset = NULL, 
                          method = c("logistic", "probit", "loglog", 
                                     "cloglog", "cauchit"), ...,
                          prior = R2(stop("'location' must be specified")), 
                          prior_counts = dirichlet(1), shape = NULL, rate = NULL, 
                          prior_PD = FALSE, 
                          algorithm = c("sampling", "meanfield", "fullrank"),
                          adapt_delta = NULL) {
  
  algorithm <- match.arg(algorithm)
  method <- match.arg(method)
  all_methods <- c("logistic", "probit", "loglog", "cloglog", "cauchit")
  link <- which(all_methods == method)
  if (!is.factor(y)) 
    stop("'y' must be a factor.")
  y_lev <- levels(y)
  J <- length(y_lev)
  y <- as.integer(y)
  if (colnames(x)[1] == "(Intercept)")
    x <- x[, -1, drop=FALSE]
  xbar <- as.array(colMeans(x))
  X <- sweep(x, 2, xbar, FUN = "-")
  cn <- colnames(X)
  decomposition <- qr(X)
  Q <- qr.Q(decomposition)
  R_inv <- qr.solve(decomposition, Q)
  X <- Q
  colnames(X) <- cn
  xbar <- c(xbar %*% R_inv)
  
  has_weights <- isTRUE(length(wt) > 0 && !all(wt == 1))
  if (!has_weights) 
    weights <- double(0)
  has_offset <- isTRUE(length(offset) > 0 && !all(offset == 0))
  if (!has_offset) 
    offset <- double(0)

  if (length(prior)) {
    regularization <- make_eta(prior$location, prior$what, K = ncol(x))
    prior_dist <- 1L
  } else {
    regularization <- 0
    prior_dist <- 0L
  }
  if (!length(prior_counts)) {
    prior_counts <- rep(1, J)
  } else {
    prior_counts <- maybe_broadcast(prior_counts$concentration, J)
  }

  if (is.null(shape)) {
    shape <- 0L
  } else {
    if (J > 2) 
      stop("'shape' must be NULL when there are more than 2 outcome categories.")
    if (!is.numeric(shape) || shape <= 0) 
      stop("'shape' must be positive")
  }
  
  if (is.null(rate)) {
    rate <- 0L
  } else {
    if (J > 2) 
      stop("'rate' must be NULL when there are more than 2 outcome categories.")
    if (!is.numeric(rate) || rate <= 0) 
      stop("'rate' must be positive")
  }

  is_skewed <- as.integer(shape > 0 & rate > 0)
  if (is_skewed && method != "logistic")
    stop("Skewed models are only supported when method = 'logistic'.")
    
  N <- nrow(X)
  K <- ncol(X)
  standata <- nlist(J, N, K, X, xbar, y, prior_PD, link, 
                    has_weights, weights, has_offset, offset,
                    prior_dist, regularization, prior_counts,
                    is_skewed, shape, rate,
                    # the rest of these are not actually used
                    has_intercept = 0L, prior_dist_for_intercept = 0L, 
                    family = 1L)
  stanfit <- stanmodels$polr
  if (J > 2) {
    pars <- c("beta", "zeta", "mean_PPD")
  } else { 
    pars <- c("zeta", "beta", if (is_skewed) "alpha", "mean_PPD")
  }
  standata$do_residuals <- isTRUE(J > 2)
  
  if (algorithm == "sampling") {
    sampling_args <- set_sampling_args(
      object = stanfit, 
      prior = prior,
      user_dots = list(...), 
      user_adapt_delta = adapt_delta, 
      data = standata, pars = pars, show_messages = FALSE)
    stanfit <- do.call(sampling, sampling_args)
  } else {
    stanfit <- rstan::vb(stanfit, pars = pars, data = standata, 
                         algorithm = algorithm, ...)
  }
  
  thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, permuted = FALSE)
  betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
  for (chain in 1:tail(dim(betas), 1)) for (param in 1:nrow(betas)) {
    stanfit@sim$samples[[chain]][[(J == 2) + param]] <- 
      if (ncol(X) > 1) betas[param, , chain] else betas[param, chain]
  }
  
  if (J > 2) {
    new_names <- c(colnames(x), 
                   paste(head(y_lev, -1), tail(y_lev, -1), sep = "|"),
                   paste("mean_PPD", y_lev, sep = ":"), 
                   "log-posterior")
  } else {
    new_names <- c("(Intercept)", 
                   colnames(x), 
                   if (is_skewed) "alpha",
                   "mean_PPD", 
                   "log-posterior")
  }
  stanfit@sim$fnames_oi <- new_names
  
  return(stanfit)
}
