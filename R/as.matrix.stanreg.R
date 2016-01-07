# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

#' Extract posterior sample
#' 
#' For models fit using MCMC (\code{algorithm="sampling"}), the posterior sample
#' is an \eqn{S} by \eqn{P} matrix, where \eqn{S} is the size of the sample --- 
#' the number of post-warmup draws from the posterior distribution of the model 
#' parameters --- and \eqn{P} is the number of parameters. If using optimization
#' (\code{"optimizing"}) or variational inference (\code{"meanfield"} or
#' \code{"fullrank"}), there is no posterior sample but rather a \eqn{1000} by
#' \eqn{P} matrix of draws from either the asymptotic multivariate Gaussian
#' sampling distribution of the parameters or the variational approximation to
#' the posterior distribution.
#' 
#' @method as.matrix stanreg
#' @export
#' 
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-regex-pars
#' @param pars An optional character vector of parameter names.
#' @param ... Ignored.
#'   
#' @return A matrix or data frame, the dimensions of which depend on \code{pars}
#'   and \code{regex_pars} (if specified), as well as the model and estimation
#'   algorithm (as described above).
#' 
#' @seealso \code{\link{stanreg-methods}}
#' 
#' @examples
#' # Extract posterior sample after MCMC
#' draws <- as.matrix(example_model)
#' 
#' # For example, we can see that the median of the draws for the intercept 
#' # is the same as the point estimate rstanarm uses
#' print(median(draws[, "(Intercept)"]))
#' print(example_model$coefficients["(Intercept)"])
#' 
#' \dontrun{
#' # Extract draws from asymptotic Gaussian sampling distribution 
#' # after optimization
#' fit <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
#' draws <- as.data.frame(fit)
#' print(colnames(draws))
#' 
#' # Extract draws from variational approximation to the posterior distribution
#' fit2 <- update(fit, algorithm = "meanfield")
#' draws <- as.data.frame(fit2)
#' print(colnames(draws))
#' }
#' 
as.matrix.stanreg <- function(x, pars = NULL, regex_pars = NULL, ...) {
  NO_DRAWS <- "No draws found."
  pars <- .collect_pars(x, pars, regex_pars)
  no_user_pars <- is.null(pars)
  if (used.optimizing(x)) {
    mat <- x$asymptotic_sampling_dist
    if (is.null(mat)) stop(NO_DRAWS, call. = FALSE)
    if (is.null(pars)) {
      dispersion <- c("sigma", "scale", "shape", "lambda", "overdispersion")
      pars <- c(names(coef(x)), # return with coefficients first
                dispersion[which(dispersion %in% colnames(mat))])
    }
  } else { # used mcmc or vb
    if (x$stanfit@mode != 0) stop(NO_DRAWS, call. = FALSE)
    posterior <- rstan::extract(x$stanfit, permuted = FALSE, inc_warmup = FALSE)
    mat <- apply(posterior, 3L, FUN = function(y) y)
    if (is.null(pars))
      pars <- grep("mean_PPD|log-posterior", colnames(mat), invert = TRUE, 
                   value = TRUE)
  }
  if (!no_user_pars) {
    badpars <- which(!pars %in% colnames(mat))
    if (length(badpars)) 
      stop("No parameter(s) ", paste(pars[badpars], collapse = ", "), 
           call. = FALSE)
  }
  mat <- mat[, pars, drop = FALSE]
  if (!is.mer(x)) return(mat) else return(unpad_reTrms(mat, columns = TRUE))
}

#' @rdname as.matrix.stanreg
#' @export
#' 
as.data.frame.stanreg <- function(x, pars = NULL, regex_pars = NULL, ...) {
  as.data.frame(as.matrix.stanreg(x, pars, regex_pars, ...))
}
