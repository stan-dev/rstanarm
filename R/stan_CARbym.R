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

#' Bayesian CAR BYM models via Stan
#'
#' Spatial regression modeling with a conditional autoregressive (CAR).
#'
#' @export
#' 

stan_CARbym <- function(formula,
                        family = c("binomial", "poisson", "gaussian"),
                        data,
                        trials = NULL,
                        W,
                        ...,
                        prior = normal(), prior_intercept = normal(),
                        prior_PD = FALSE,
                        algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                        adapt_delta = NULL,
                        QR = FALSE) {
  stan_function <- "stan_CARbym"
  if (!requireNamespace("CARBayes", quietly = TRUE))
    stop("Please install the CARBayes package before using 'stan_CARbym.")
  mc <- match.call(expand.dots = FALSE)
  algorithm <- match.arg(algorithm)
  family <- match.arg(family)
  mf <- model.frame(mc, data)
  Y <- array1D_check(model.response(mf, type = "any"))
  X <- model.matrix(formula, data)
  stanfit <- stan_spatial.fit(x = X, y = Y, w = W,
                              trials = trials,
                              family = family,
                              stan_function,
                              ...,
                              prior = prior,
                              prior_intercept = prior_intercept,
                              prior_PD = prior_PD,
                              algorithm = algorithm, adapt_delta = adapt_delta, 
                              QR = QR)
  if (family == "binomial")
    family <- binomial(link = "logit")
  else if (family == "poisson")
    family <- poisson(link = "log")
  else if (family == "gaussian")
    family <- gaussian(link = "identity")
  fit <- nlist(stanfit,
               algorithm,
               data,
               x = X, y = Y,
               family,
               formula,
               model = mf, 
               call = match.call(),
               stan_function = "stan_CARbym")
  out <- stanreg(fit)
  browser()
  structure(out, class = c("stanreg", "car"))
}
