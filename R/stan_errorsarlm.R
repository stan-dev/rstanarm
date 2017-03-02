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

#' Bayesian spatial simultaneous autoregressive error model estimation via Stan
#'
#' Note this doc is incomplete!
#'
#' @export
#' 
#' @param listw Spatial weights as a "listw" object. This can be constructed from
#' a variety of formats using the appropriate functions in the \code{spdep}
#' package (e.g. \code{mat2listw} transforms a "matrix" class object to a
#' "listw" class object).
#' @param  prior_rho Prior on spatial autocorrelation term.
#' @param prior_intercept Prior on intercept of linear predictor.
#' 
#' @examples 
#' ### Spatial AR error Simulation
#' N <- 10
#' W_bin <- matrix(rep(0, N * N), nrow = N)
#' W_bin[lower.tri(W_bin)] <- rbinom(choose(N,2), 1, 0.5)
#' W_bin <- W_bin + t(W_bin)
#' W <- apply(W_bin, 2, function(x){x/rowSums(W_bin)})
#' I <- diag(N)
#' lambda <- 0.5
#' sigma <- 0.3
#' Sigma <- solve((I - lambda * W) %*% (I - lambda * t(W))) * sigma
#' X <- cbind(rep(1,N),rnorm(N, 0, 1), rnorm(N, 3, 1))
#' beta <- c(3, 2.5, -1.5)
#' mu <- X %*% beta
#' y <- c(mvtnorm::rmvnorm(1, mu, Sigma))
#' lw <- spdep::mat2listw(W)
#' 
#' dat <- data.frame(cbind(y, X[,-1]))
#' names(dat) <- c("y","x1","x2")
#' 
#' fit <- stan_errorsarlm(y ~ x1 + x2, data = dat, listw = lw, cores = 4)

stan_errorsarlm <- function(formula, data, listw, type = "lag", ...,
                          prior_rho = beta(), prior_intercept = NULL,
                          algorithm = c("sampling", "optimizing", "meanfield", "fullrank"), 
                          adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  sp_model <- "errorsarlm"
  
  if (!requireNamespace("spdep", quietly = TRUE))
    stop("Please install the spdep package before using 'stan_lagsarlm'.")
  algorithm <- match.arg(algorithm)
  validate_glm_formula(formula)
  
  mc <- match.call(expand.dots = FALSE)
  
  # NULLify any Stan specific arguments in mc
  mc$prior_rho <- mc$prior_intercept <- mc$algorithm <- mc$adapt_delta <- 
    mc$QR <- mc$sparse <- NULL
  
  # mc$drop.unused.levels <- TRUE  # drop stuff in ... so quote evaluates
  mc[[1L]] <- quote(spdep::lagsarlm)
  mc$... <- NULL
  sp <- suppressWarnings(eval(mc, parent.frame()))
  
  Y <- array1D_check(sp$y)
  X <- sp$X
  W <- spdep::listw2mat(listw)
  
  stanfit <- stan_sp.fit(y = Y, x = X, w = W, ..., sp_model = sp_model,
                         prior_rho = prior_rho, prior_intercept = prior_intercept, 
                         algorithm = algorithm, adapt_delta = adapt_delta, QR = QR, sparse = sparse)
  
  
  fit <- nlist(stanfit, algorithm, data, family = gaussian(),
               x = X, y = Y, formula, call = match.call(), sp_model)
  
  
  out <- stanreg(fit)
  
  class(out) <- c("stanreg", "spatial")
  
  return(out)
}