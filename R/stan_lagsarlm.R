# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
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

#' Bayesian spatial simultaneous autogregressive lag model estimation
#'
#' @export

stan_lagsarlm <- function(formula, data, listw, type = "lag", ...,
                          prior_rho = NULL, prior_intercept = NULL,
                          algorithm = c("sampling", "optimizing", "meanfield", "fullrank"), 
                          adapt_delta = NULL) {
  sp_model <- "lagsarlm"
  
  if (!requireNamespace("spdep", quietly = TRUE))
    stop("Please install the spdep package before using 'stan_lagsarlm'.")
  algorithm <- match.arg(algorithm)
  validate_glm_formula(formula)

  mc <- match.call(expand.dots = FALSE)
  
  # NULLify any Stan specific arguments in mc
  mc$prior_rho <- mc$prior_intercept <- mc$algorithm <- mc$adapt_delta <- NULL
  
  # mc$drop.unused.levels <- TRUE  # drop stuff in ... so quote evaluates
  mc[[1L]] <- quote(spdep::lagsarlm)
  mc$... <- NULL
  sp <- suppressWarnings(eval(mc, parent.frame()))
  
  Y <- array1D_check(sp$y)
  X <- sp$X
  W <- spdep::listw2mat(listw)
  
  stanfit <- stan_sp.fit(y = Y, x = X, w = W, ..., sp_model = sp_model,
                          prior_rho = prior_rho, prior_intercept = prior_intercept, 
                          algorithm = algorithm, adapt_delta = adapt_delta)
  
  
  fit <- nlist(stanfit, algorithm, data, family = gaussian(),
               x = X, y = Y, formula, call = match.call())
  
  
  out <- stanreg(fit)

  structure(out, class = c("stanreg", "sarlm"))

  return(out)
}