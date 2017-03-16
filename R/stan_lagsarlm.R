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

#' Bayesian spatial simultaneous autoregressive lag model estimation via Stan
#'
#' Spatial autoregressive (SAR) models with spatial dependence modeled through
#' a spatially lagged dependent variable. The model takes the form:
#' \deqn{y = \rho W y + X \beta + e}
#' where \eqn{\rho} is the spatial autoregressive lag coefficient and \eqn{W}
#' is a spatial weight matrix.
#' 
#' @export
#' @templateVar fun stan_lagsarlm
#' @templateVar fitfun stan_sp.fit
#' @template args-x-y
#' @param listw Spatial weights as a "listw" object. This can be constructed from
#' a variety of formats using the appropriate functions in the \code{spdep}
#' package (e.g. \code{spdep::mat2listw} transforms a "matrix" class object to a
#' "listw" class object).
#' @param  prior_aux Prior on spatial autocorrelation term.
#'   \code{prior_aux} can be set to \code{beta}.
#'   See the \link[=priors]{priors help page} for details.
#'   To omit a prior on the intercept ---i.e., to use a flat
#'   (improper) uniform prior--- \code{prior_aux} can be set to
#'   \code{NULL}.
#' @param prior_intercept The prior distribution for the intercept. 
#'   \code{prior_intercept} can be a call to \code{normal}, \code{student_t} or 
#'   \code{cauchy}. See the \link[=priors]{priors help page} for details on 
#'   these functions. To omit a prior on the intercept ---i.e., to use a flat
#'   (improper) uniform prior--- \code{prior_intercept} can be set to
#'   \code{NULL}.
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-dots
#' @template return-stanreg-object
#' @template return-stanfit-object
#' 
#' @details The \code{stan_lagsarlm} function is similar in syntax to 
#'   \code{\link[spdep]{lagsarlm}} but rather than performing maximum likelihood 
#'   estimation, full Bayesian estimation is performed (if \code{algorithm} is 
#'   \code{"sampling"}) via MCMC. The Bayesian model adds priors on the intercept
#'   and the spatial autocorrelation coefficient. The \code{stan_lagsarlm}
#'   function calls the workhorse \code{stan_sp.fit} function.
#' 
#' @examples 
#' \donttest{
#' ### Spatial AR lag Simulation
#' path <- system.file(package = "rstanarm", "data/spatial")
#' sim_grid <- rgdal::readOGR(path, layer = "grid_map")
#' N <- nrow(as.data.frame(sim_grid))
#' I <- diag(N)
#' lambda <- 0.8
#' sigma <- 0.3
#' W <- spdep::nb2mat(spdep::poly2nb(sim_grid, queen = TRUE), style = "W", zero.policy = TRUE)
#' Sigma <- solve(I - lambda * W) %*% t(solve(I - lambda * W)) * sigma
#' X <- cbind(rep(1,N),rnorm(N, 0, 1), rnorm(N, 3, 1))
#' beta <- c(3, 2.5, 1.5)
#' mu <- solve(I - lambda * W) %*% X %*% beta
#' # y <- c(mvtnorm::rmvnorm(1, mu, Sigma))
#' y <- rstanarm:::rmultinorm(1, mu, Sigma)
#' lw <- spdep::mat2listw(W)
#' dat <- data.frame(cbind(y, X[,-1]))
#' names(dat) <- c("y","x1","x2")
#' fit <- stan_lagsarlm(y ~ x1 + x2, data = dat, listw = lw, cores = 4, iter = 100)
#' print(fit, digits = 2)
#' }

stan_lagsarlm <- function(formula, data, listw, type = "lag", ...,
                          prior_aux = beta(), prior_intercept = normal(),
                          algorithm = c("sampling", "optimizing", "meanfield", "fullrank"), 
                          adapt_delta = NULL, QR = TRUE, sparse = FALSE) {
  sp_model <- "lagsarlm"

  if (!requireNamespace("spdep", quietly = TRUE))
    stop("Please install the spdep package before using 'stan_lagsarlm'.")
  if (!("listw" %in% class(listw)))
    stop("Spatial weights must be a listw object. See the spdep package documentation for details.")
  
  algorithm <- match.arg(algorithm)
  validate_glm_formula(formula)
  mc <- match.call(expand.dots = FALSE)
  # NULLify any Stan specific arguments in mc
  mc$prior_aux <- mc$prior_intercept <- mc$algorithm <- mc$adapt_delta <- 
    mc$QR <- mc$sparse <- NULL
  
  # mc$drop.unused.levels <- TRUE  # drop stuff in ... so quote evaluates
  mc[[1L]] <- quote(spdep::lagsarlm)
  mc$... <- NULL
  sp <- suppressWarnings(eval(mc, parent.frame()))

  Y <- array1D_check(sp$y)
  X <- sp$X
  W <- spdep::listw2mat(listw)
  
  if (!(all(rowSums(W) == rep(1,nrow(W))))) {
    warning("Spatial weights are not row normalized.", call. = FALSE)
  }

  stanfit <- stan_sp.fit(y = Y, x = X, w = W, ..., sp_model = sp_model,
                          prior_aux = prior_aux, prior_intercept = prior_intercept, 
                          algorithm = algorithm, adapt_delta = adapt_delta, QR = QR, sparse = sparse)
  
  
  fit <- nlist(stanfit, algorithm, data, family = gaussian(),
               x = X, y = Y, formula, call = match.call(), sp_model)
  
  
  out <- stanreg(fit)

  class(out) <- c("stanreg", "spatial")
  
  return(out)
}