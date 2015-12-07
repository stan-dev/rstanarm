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
                         prior_intercept = NULL, prior_PD = FALSE, 
                         algorithm = c("sampling", "meanfield", "fullrank"),
                         adapt_delta = NULL) {
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
  xbar <- colMeans(x)
  x <- sweep(x, 2, xbar, FUN = "-")
  b <- coef(ols)
  b[is.na(b)] <- 0.0
  cn <- colnames(x)
  decomposition <- qr(x)
  Q <- qr.Q(decomposition)
  R <- qr.R(decomposition)
  R_inv <- qr.solve(decomposition, Q)
  x <- Q
  colnames(x) <- cn
  xbarR_inv <- array(c(xbar %*% R_inv), c(J,K))
  if (has_intercept == 1) Rb <- array(R %*% b[-1], c(J,K))
  else Rb <- array(R %*% b, c(J,K))

  SSR <- array(crossprod(residuals(ols))[1], J)
  s_Y <- array(sd(y), J)
  if (isTRUE(all.equal(matrix(0, J, K), xbar))) center_y <- mean(y)
  else center_y <- 0
  ybar <- array(mean(y), J)

  if (length(prior) == 0) {
    prior_dist <- 0L
    eta <- 0
  }
  else {
    prior_dist <- 1L
    eta <- prior$eta <- make_eta(prior$location, prior$what, K = K)
  }
  if (length(prior_intercept) == 0) {
    prior_dist_for_intercept <- 0L
    prior_mean_for_intercept <- 0
    prior_scale_for_intercept <- 0
  }
  else {
    if (!identical(prior_intercept$dist, "normal"))
      stop("'prior_intercept' must be 'NULL' or a call to 'normal'")
    prior_dist_for_intercept <- 1L
    prior_mean_for_intercept <- prior_intercept$location
    prior_scale_for_intercept <- prior_intercept$scale
    if (is.null(prior_scale_for_intercept)) {
      prior_scale_for_intercept <- 0
    }
  }
  dim(R_inv) <- c(J, dim(R_inv))
  
  # initial values
  R2 <- array(1 - SSR[1] / ((N - 1) * var(y)), J)
  log_omega <- array(0, ifelse(prior_PD == 0, J, 0))
  init_fun <- function(chain_id) {
    out <- list(R2 = R2, log_omega = log_omega)
    return(out)
  }
  algorithm <- match.arg(algorithm)
  stanfit <- stanmodels$lm
  standata <- nlist(K, has_intercept, prior_dist,
                    prior_dist_for_intercept, 
                    prior_mean_for_intercept,
                    prior_scale_for_intercept,
                    prior_PD, eta, J, N, xbarR_inv,
                    ybar, center_y, s_Y, Rb, SSR, R_inv)
  pars <- c(if (has_intercept) "alpha", "beta", "sigma", 
            if (prior_PD == 0) "log_omega", "R2", "mean_PPD")
  if (algorithm == "optimizing") {
    stop("'optimizing' is not a supported estimation technique for this model")
    opt <- optimizing(stanfit, data = standata, init = init_fun,
                      draws = 1000, ...) # yields corner solution
    new_names <- names(opt$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    new_names[mark] <- cn
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    new_names[new_names == "sigma[1]"] <- "sigma"
    new_names[new_names == "log_omega[1]"] <- "log-fit_ratio"
    new_names[new_names == "R2[1]"] <- "R2"
    names(opt$par) <- new_names
    colnames(opt$theta_tilde) <- new_names
    opt$stanfit <- suppressMessages(sampling(stanfit, data = standata, chains = 0))
    return(opt)
  }
  else if (algorithm %in% c("meanfield", "fullrank")) {
    stanfit <- rstan::vb(stanfit, data = standata, pars = pars,
                         algorithm = algorithm, ...)
  }
  else {
    sampling_args <- set_sampling_args(
      object = stanfit, 
      prior = prior,
      user_dots = list(...), 
      user_adapt_delta = adapt_delta, 
      init = init_fun, data = standata, pars = pars, show_messages = FALSE)
    stanfit <- do.call(sampling, sampling_args)
  }
  new_names <- c(if (has_intercept) "(Intercept)", colnames(x), "sigma", 
                 if (prior_PD == 0) "log-fit_ratio", "R2", "mean_PPD", "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  return(stanfit)
}

#' @rdname stan_lm
#' @export
stan_lm.fit <- function(x, y, offset = NULL, singular.ok = TRUE, ...,
                        prior = R2(stop("'location' must be specified")), 
                        prior_intercept = NULL, prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta = NULL) {

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf[[1L]] <- as.name("stan_lm.wfit")
  mf$w <- as.name("NULL")
  return(eval(mf, parent.frame()))
}
