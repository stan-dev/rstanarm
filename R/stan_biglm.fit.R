# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016, 2017 Trustees of Columbia University
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
#' @param b A numeric vector of OLS coefficients, excluding the intercept
#' @param R A square upper-triangular matrix from the QR decomposition of the 
#'   design matrix, excluding the intercept
#' @param SSR A numeric scalar indicating the sum-of-squared residuals for OLS
#' @param N A integer scalar indicating the number of included observations
#' @param has_intercept A logical scalar indicating whether to add an intercept 
#'   to the model when estimating it.
#' @param importance_resampling Logical scalar indicating whether to use 
#'   importance resampling when approximating the posterior distribution with
#'   a multivariate normal around the posterior mode, which only applies
#'   when \code{algorithm} is \code{"optimizing"} but defaults to \code{TRUE}
#'   in that case
#' @param keep_every Positive integer, which defaults to 1, but can be higher
#'   in order to thin the importance sampling realizations and also only
#'   apples when \code{algorithm} is \code{"optimizing"} but defaults to
#'   \code{TRUE} in that case
#' @examples
#' # create inputs
#' ols <- lm(mpg ~ wt + qsec + am, data = mtcars, # all row are complete so ...
#'           na.action = na.exclude)              # not necessary in this case
#' b <- coef(ols)[-1]
#' R <- qr.R(ols$qr)[-1,-1]
#' SSR <- crossprod(ols$residuals)[1]
#' not_NA <- !is.na(fitted(ols))
#' N <- sum(not_NA)
#' xbar <- colMeans(mtcars[not_NA,c("wt", "qsec", "am")])
#' y <- mtcars$mpg[not_NA]
#' ybar <- mean(y)
#' s_y <- sd(y)
#' post <- stan_biglm.fit(b, R, SSR, N, xbar, ybar, s_y, prior = R2(.75),
#'                        # the next line is only to make the example go fast
#'                        chains = 1, iter = 500, seed = 12345)
#' cbind(lm = b, stan_lm = rstan::get_posterior_mean(post)[13:15,]) # shrunk
#' 
stan_biglm.fit <- function(b, R, SSR, N, xbar, ybar, s_y, has_intercept = TRUE, ...,
                           prior = R2(stop("'location' must be specified")), 
                           prior_intercept = NULL, prior_PD = FALSE, 
                           algorithm = c("sampling", "meanfield", "fullrank", "optimizing"),
                           adapt_delta = NULL,
                           importance_resampling = TRUE,
                           keep_every = 1) {
  
  if (prior_PD && is.null(prior_intercept)) {
    msg <- "The default flat prior on the intercept is not recommended when 'prior_PD' is TRUE."
    warning(msg, call. = FALSE, immediate. = TRUE)
    warning(msg, call. = FALSE, immediate. = FALSE)
  }
  
  J <- 1L
  N <- array(N, c(J))
  K <- ncol(R)
  cn <- names(xbar)
  if (is.null(cn)) cn <- names(b)
  R_inv <- backsolve(R, diag(K))
  JK <- c(J, K)
  xbarR_inv <- array(c(xbar %*% R_inv), JK)
  Rb <- array(R %*% b, JK)
  SSR <- array(SSR, J)
  s_Y <- array(s_y, J)
  center_y <- if (isTRUE(all.equal(matrix(0, J, K), xbar))) ybar else 0.0
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
    
    # also add scale back to prior_intercept to pass to summarize_lm_prior later
    prior_intercept$scale <- prior_scale_for_intercept
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
  pars <- c(if (has_intercept) "alpha", "beta", "sigma", 
            if (prior_PD == 0) "log_omega", "R2", "mean_PPD")
  algorithm <- match.arg(algorithm)
  if (algorithm == "optimizing") {
    optimizing_args <- list(...)
    if (is.null(optimizing_args$draws)) optimizing_args$draws <- 1000L
    optimizing_args$object <- stanfit
    optimizing_args$data <- standata
    optimizing_args$constrained <- TRUE
    optimizing_args$importance_resampling <- importance_resampling
    if (is.null(optimizing_args$tol_rel_grad)) 
      optimizing_args$tol_rel_grad <- 10000L
    out <- do.call(optimizing, args = optimizing_args)
    check <- check_stanfit(out)
    if (!isTRUE(check)) return(standata)
    if (K == 1)
        out$theta_tilde[,'R2[1]'] <- (out$theta_tilde[,'R2[1]']) ^ 2
    pars_idx <- unlist(sapply(1:length(pars), function(i) {
      which(grepl(paste('^', pars[i], sep=''), names(out$par)))
    }))
    nrows <- dim(out$theta_tilde)[1]
    out$theta_tilde <- out$theta_tilde[,pars_idx]
    dim(out$theta_tilde) <- c(nrows, length(pars_idx))
    new_names <- c(if (has_intercept) "(Intercept)", cn, "sigma", 
                   if (prior_PD == 0) "log-fit_ratio", 
                   "R2", "mean_PPD")
    colnames(out$theta_tilde) <- new_names
    if (optimizing_args$draws > 0) { # begin: psis diagnostics and importance resampling
        lr <- out$log_p-out$log_g
        lr[lr == -Inf] <- -800
        p <- suppressWarnings(loo::psis(lr, r_eff = 1))
        p$log_weights <- p$log_weights - log_sum_exp(p$log_weights)
        theta_pareto_k <- suppressWarnings(apply(out$theta_tilde, 2L, function(col) {
          if (all(is.finite(col))) loo::psis(log1p(col ^ 2) / 2 + lr, r_eff = 1)$diagnostics$pareto_k else NaN
        }))
        ## todo: change fixed threshold to an option
        if (any(theta_pareto_k > 0.7, na.rm = TRUE)) {
            warning("Some Pareto k diagnostic values are too high. Resampling disabled.",
                    "Decreasing tol_rel_grad may help if optimization has terminated prematurely.", 
                    " Otherwise consider using sampling instead of optimizing.", call. = FALSE, immediate. = TRUE)
            importance_resampling <- FALSE
        } else if (any(theta_pareto_k > 0.5, na.rm = TRUE)) { 
            warning("Some Pareto k diagnostic values are slightly high.",
                    " Increasing the number of draws or decreasing tol_rel_grad may help.", 
                    call. = FALSE, immediate. = TRUE)
        }
        out$psis <- nlist(pareto_k = p$diagnostics$pareto_k, n_eff = p$diagnostics$n_eff / keep_every)
    } else {
      theta_pareto_k <- rep(NaN, length(new_names))
      importance_resampling <- FALSE
    }
    if (importance_resampling) {  
      ir_idx <- .sample_indices(exp(p$log_weights), 
                                n_draws = ceiling(optimizing_args$draws / keep_every))
      out$theta_tilde <- out$theta_tilde[ir_idx,]
      out$ir_idx <- ir_idx
      ## SIR mcse and n_eff
      w_sir <- as.numeric(table(ir_idx)) / length(ir_idx)
      mcse <- apply(out$theta_tilde[!duplicated(ir_idx),], 2L, function(col) {
        if (all(is.finite(col))) sqrt(sum(w_sir ^ 2 * (col-mean(col)) ^ 2)) 
        else NaN
      })
      n_eff <- round(apply(out$theta_tilde[!duplicated(ir_idx),], 2L, var) / (mcse^2), digits = 0)
    } else {
      out$ir_idx <- NULL
      mcse <- rep(NaN, length(theta_pareto_k))
      n_eff <- rep(NaN, length(theta_pareto_k))
    }
    out$diagnostics <- cbind(mcse, theta_pareto_k, n_eff)
    colnames(out$diagnostics) <- c("mcse", "khat", "n_eff")
    ## end: psis diagnostics and SIR
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, 
                                             chains = 0))
    prior_info <- summarize_lm_prior(prior, prior_intercept)
    return(structure(out, prior.info = prior_info))
  } else if (algorithm %in% c("meanfield", "fullrank")) {
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
  check <- check_stanfit(stanfit)
  if (!isTRUE(check)) return(standata)
  if (K == 1)
    stanfit@sim$samples <- lapply(stanfit@sim$samples, FUN = function(x) {
      x$`R2[1]` <- (x$`R2[1]`)^2
      return(x)
    })
  new_names <- c(if (has_intercept) "(Intercept)", cn, "sigma", 
                 if (prior_PD == 0) "log-fit_ratio", 
                 "R2", "mean_PPD", "log-posterior")
  stanfit@sim$fnames_oi <- new_names

  prior_info <- summarize_lm_prior(prior, prior_intercept)
  structure(stanfit, prior.info = prior_info)
}


# internal ----------------------------------------------------------------

# Create "prior.info" attribute needed for prior_summary()
#
# @param prior, prior_intercept User's prior and prior_intercept specifications
# @return A named list with elements 'prior' and 'prior_intercept' containing 
#   the values needed for prior_summary
summarize_lm_prior <- function(prior, prior_intercept) {
  flat <- !length(prior)
  flat_int <- !length(prior_intercept)
  
  list(
    prior = list(
      dist = ifelse(flat, NA, "R2"),
      location = ifelse(flat, NA, prior$location),
      what = ifelse(flat, NA, prior$what)
    ), 
    prior_intercept = list(
      dist = ifelse(flat_int, NA, "normal"),
      location = ifelse(flat_int, NA, prior_intercept$location),
      scale = ifelse(flat_int, NA, prior_intercept$scale)
    )
  )
}
