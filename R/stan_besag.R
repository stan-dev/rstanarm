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

#' Bayesian CAR Intrinsic Autoregressive models via Stan
#'
#' Spatial regression modeling with an intrinsic conditional autoregressive (ICAR) prior.

#' @rdname stan_besag
#' @export
#' 
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template see-also
#' @template args-formula-data-subset
#' @template args-same-as
#' @template args-same-as-rarely
#' @template args-x-y
#' @template args-dots
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' 
#' @param family Distribution associated with the outcome. Gaussian, Binomial,
#' and Poisson families are supported.
#' @param trials If \code{family = binomial()} then a vector of trials (equal
#' in length to the outcome) must be declared.
#' @param W An N-by-N spatial weight matrix.
#' @param prior_rho The prior distribution on the variance of the non-centered structured (spatial) effect.
#' @param prior_sigma The prior distribution on the standard deviation of the outcome if \code{family = gaussian()} is declared.
#' 
#' @details The \code{stan_besag} model is similar to the analogous model in R-INLA. However, instead of using the integrated Laplace approximation (INLA) method, full Bayesian estimation is performed (if \code{algorithm} is \code{"sampling"}) via MCMC. The model includes priors on the intercept, regression coefficients, and the relevant scale parameters. The \code{stan_besag} function calls the workhorse \code{stan_spatial.fit} function, but it is also possible to call the latter directly.
#'   
#' @seealso The vignette for \code{stan_besag}.
#' 
#' @references Riebler, A., Sorbye, S.H., Simpson, D., Rue, H. (2016). An intuitive Bayesian spatial model for disease mapping that accounts for scaling. arXiv preprint	arXiv:1601.01180.
#' 
#' @examples 
#' ### Simulated Data on a Lattice
#' 
#' data("lattice10", package = "rstanarm")
#' 
#' # plot GMRF
#' var_range_gmrf <- seq(min(grid_sim@data$gmrf), max(grid_sim@data$gmrf), length = 50)
#' spplot(grid_sim, "gmrf", at = var_range_gmrf, main = expression(paste(phi, " (GMRF)")),
#'        col.regions = colorRampPalette(c("#ef8a62", "#f7f7f7", "#67a9cf"))(50))
#' 
#' # Convert a spatial polygon to an N-by-N weight matrix
#' sp2weightmatrix <- function(spatialpolygon) {
#'   spdep::nb2mat(spdep::poly2nb(spatialpolygon, queen = TRUE), style = "B", zero.policy = TRUE)
#' }
#' 
#' # convert spatial object to neighborhood matrix
#' W <- sp2weightmatrix(grid_sim)
#' # W_sparse <- Matrix(W, sparse = TRUE, dimnames = list(NULL,NULL)) 
#' 
#' # simulate predictor/outcome
#' x <- rnorm(nrow(W), 3, 1)
#' spatial_data <- data.frame(x, phi = grid_sim@data$gmrf)
#' spatial_data$y_gauss <- rnorm(nrow(W), 0 + 0.4 * x + spatial_data$phi, 1)
#' 
#' # fit the model
#' fit_besag <- stan_besag(y_gauss ~ 1 + x + I(x^2), data = spatial_data, W = W, iter = 300, chains = 4)
#' pp_check(fit_besag)
#'

stan_besag <- function(formula,
                     family = gaussian(),
                     data,
                     trials = NULL,
                     W,
                     ...,
                     prior = normal(), prior_intercept = normal(),
                     prior_sigma = NULL, prior_rho = normal(),
                     prior_PD = FALSE,
                     algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                     adapt_delta = NULL,
                     QR = FALSE) {
  stan_function <- "stan_besag"
  if (!requireNamespace("INLA", quietly = TRUE))
    stop(paste("Please install and load the INLA package before using", stan_function))
  mc <- match.call(expand.dots = FALSE)
  algorithm <- match.arg(algorithm)
  family <- validate_family(family)
  mf <- model.frame(mc, data)
  Y <- array1D_check(model.response(mf, type = "any"))
  X <- model.matrix(formula, data)
  
  stanfit <- stan_spatial.fit(x = X, y = Y, w = W,
                              trials = trials,
                              family = family,
                              stan_function = stan_function,
                              ...,
                              prior = prior,
                              prior_intercept = prior_intercept,
                              prior_sigma = prior_sigma,
                              prior_rho = prior_rho,
                              prior_PD = prior_PD,
                              algorithm = algorithm, adapt_delta = adapt_delta, 
                              QR = QR)
  fit <- nlist(stanfit,
               algorithm,
               data,
               x = X, y = Y,
               family,
               formula,
               model = mf, 
               call = match.call(),
               stan_function = stan_function)
 
  if (family$family == "binomial") {
    fit$family <- binomial(link = "logit")
    fit$trials <- trials
  }
  out <- stanreg(fit)
  structure(out, class = c("stanreg", "car"))
}
