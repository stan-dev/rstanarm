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

#' Bayesian spatial CAR BYM variant models via Stan
#'
#' Spatial regression modeling with a variant of the Besag, York, Mollie (BYM)
#' conditional autoregressive (CAR) prior that accounts for scaling.
#' 
#' @export
#' 
#' @templateVar fun stan_bym2
#' @templateVar fitfun stan_spatial.fit
#' @templateVar pkg rstanarm
#' @templateVar pkgfun stan_glm
#' @templateVar sameargs family 
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template args-formula-data-subset
#' @template args-same-as
#' @template args-x-y
#' @template args-dots
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' 
#' @param trials If \code{family = binomial()} then a vector of trials (equal in
#'   length to the outcome) must be declared.
#' @param W An N-by-N spatial weight matrix.
#' @param prior_structured The prior on the marginal variance contribution of the
#'   structured (spatial) and unstructured (random) effect.
#' @param prior_mixing The prior on the proportion of the marginal variance that is
#'   explained by the structured (spatial) effect. The hyperparameter \code{rho}
#'   is on the unit interval so users have the option of declaring a Beta prior 
#'   distribution or a flat prior. A prior distribution with most of the mass
#'   around 1 is analogous to the prior belief that there exists a strong
#'   spatial relationship on the graph.
#' @param order Order of the spatial random walk. Specifying \code{order = 2}
#'   will smooth the spatial variation. The default is \code{order = 1}.
#'   
#' @details The \code{stan_bym2} model is similar to the BYM2 model in R-INLA.
#'   However, instead of using the integrated nested Laplace approximation 
#'   (INLA) method, full Bayesian estimation is performed (if \code{algorithm} 
#'   is \code{"sampling"}) via MCMC. The model includes priors on the intercept,
#'   regression coefficients, spatial mixing parameter, overall spatial
#'   variation, and any applicable auxiliary parameters. The \code{stan_bym2}
#'   function calls the workhorse \code{stan_spatial.fit} function, but it is
#'   also possible to call the latter directly.
#' 
#' @seealso The vignette for \code{stan_bym2}.
#' 
#' @references Riebler, A., Sorbye, S.H., Simpson, D., Rue, H. (2016). An
#'   intuitive Bayesian spatial model for disease mapping that accounts for
#'   scaling. arXiv preprint	arXiv:1601.01180.
#'   
#'   Besag, J., York, J. and Mollié, A. (1991). Bayesian image restoration, with
#'   two applications in spatial statistics. Annals of the Institute of
#'   Statistical Mathematics. Vol. 43, No. 01, p1-20.
#'   
#'   Simpson, D., Rue, H., Martins, T.G., Riebler, A. and Sorbye, S.H. (2015).
#'   Penalising model component complexity: A principled, practical approach to
#'   constructing priors. arXiv preprint	arXiv:1403.4630.
#' 
#' @examples 
#' \dontrun{
#' ### Simulated Data on a Lattice
#' 
#' data("lattice", package = "rstanarm")
#' 
#' # plot GMRF
#' grid_sim <- grid_sim15
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
#' phi <- grid_sim@data$gmrf
#' theta <- rnorm(nrow(W), 0, 1)
#' tau <- 1
#' rho <- 0.7
#' psi <- (1/tau) * (sqrt(1-rho)*theta + sqrt(rho)*phi)
#' y <- rnorm(nrow(W), 0 + 0.4 * x + psi, 1)
#' 
#' # fit the model
#' fit_bym2 <- stan_bym2(y ~ 1 + x, data = data.frame(y=y,x=x),
#'                       W = W, iter = 1e3, chains = 4, cores = 2)
#' fit_bym2
#' pp_check(fit_bym2)
#' }
#' 

stan_bym2 <- function(formula,
                        family = gaussian(),
                        data,
                        trials = NULL,
                        W,
                        order = 1,
                        ...,
                        prior = normal(), prior_intercept = normal(),
                        prior_structured = normal(), prior_mixing = beta(0.5,0.5), prior_aux = NULL,
                        prior_PD = FALSE,
                        algorithm = c("sampling", "meanfield", "fullrank"),
                        adapt_delta = NULL,
                        QR = FALSE) {
  stan_function <- "stan_bym2"
  mc <- match.call(expand.dots = FALSE)
  algorithm <- match.arg(algorithm)
  family <- validate_family(family)
  mf <- model.frame(mc, data)
  mt <- terms(formula, data = data)
  Y <- array1D_check(model.response(mf, type = "any"))
  X <- model.matrix(formula, data)
  
  if (!is.null(prior_mixing)) {
    if (prior_mixing$dist != "beta")
      stop("'prior_mixing' must be either beta() or NULL.") 
  }
  
  stanfit <- stan_spatial.fit(x = X, y = Y, w = W,
                              trials = trials,
                              family = family,
                              stan_function = stan_function,
                              order = order,
                              ...,
                              prior = prior,
                              prior_intercept = prior_intercept,
                              prior_aux = prior_aux, prior_rho = prior_mixing, prior_tau = prior_structured,
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
               terms = mt,
               call = match.call(),
               stan_function = stan_function)

  if (family$family == "binomial") {
    fit$family <- binomial(link = "logit")
    fit$trials <- trials
  }
  out <- stanreg(fit)
  structure(out, class = c("stanreg", "car"))
}
