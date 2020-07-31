# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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
  if (NCOL(y) > 1) {
    stop("Multivariate responses not supported yet.")
  }
  
  if (colnames(x)[1L] == "(Intercept)") {
    has_intercept <- 1L
    x <- x[, -1L, drop = FALSE]
    if (NCOL(x) == 0L) {
      stop("'stan_lm' is not suitable for estimating a mean.",
           "\nUse 'stan_glm' with 'family = gaussian()' instead.", 
           call. = FALSE)
    }
  } else {
    has_intercept <- 0L
  }
  
  if (nrow(x) < ncol(x)) {
    stop("stan_lm with more predictors than data points is not yet enabled.", 
         call. = FALSE)
  }
  
  # allow prior_PD even if no y variable
  if (is.null(y)) {
    if (!prior_PD) {
      stop("Outcome variable must be specified if 'prior_PD' is not TRUE.")
    } else {
      y <- fake_y_for_prior_PD(N = NROW(x), family = gaussian())
    }
  }
  
  xbar <- colMeans(x)
  x <- sweep(x, 2L, xbar, FUN = "-")
  ybar <- mean(y)
  y <- y - ybar
  ols <- if (length(w) == 0) lm.fit(x, y) else lm.wfit(x, y, w)
  b <- coef(ols)
  NAs <- is.na(b)
  if (any(NAs) && singular.ok) {
    x <- x[,!NAs, drop = FALSE]
    xbar <- xbar[!NAs]
    ols <- lsfit(x, y, w, intercept = FALSE)
    b <- coef(ols)
  } else {
    b[NAs] <- 0.0
  }
  
  if (!is.null(w)) {
    x <- sqrt(w) * x
  }

  return(stan_biglm.fit(b, R = qr.R(ols$qr), SSR = crossprod(residuals(ols))[1], 
                        N = nrow(x), xbar = xbar, ybar = ybar, s_y = sd(y),
                        has_intercept = has_intercept, ...,
                        prior = prior, prior_intercept = prior_intercept,
                        prior_PD = prior_PD, algorithm = algorithm, 
                        adapt_delta = adapt_delta))
  
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
