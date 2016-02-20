# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-pars
#' @template args-regex-pars
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
#' print(example_model$coefficients[["(Intercept)"]])
#' 
#' \dontrun{
#' # Extract draws from asymptotic Gaussian sampling distribution 
#' # after optimization
#' fit <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
#' draws <- as.data.frame(fit)
#' print(colnames(draws))
#' print(nrow(draws)) # 1000 draws are taken
#' 
#' # Extract draws from variational approximation to the posterior distribution
#' fit2 <- update(fit, algorithm = "meanfield")
#' draws <- as.data.frame(fit2, pars = "wt")
#' print(colnames(draws))
#' print(nrow(draws)) # 1000 draws are taken
#' }
#' 
as.matrix.stanreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  STOP_no_draws <- function() stop("No draws found.", call. = FALSE)
  pars <- collect_pars(x, pars, regex_pars)
  user_pars <- !is.null(pars)
  
  if (used.optimizing(x)) {
    mat <- x$asymptotic_sampling_dist
    if (is.null(mat)) 
      STOP_no_draws()
    if (!user_pars) {
      dispersion <- c("sigma", "scale", "shape", "lambda", "overdispersion")
      pars <- c(names(coef(x)), # return with coefficients first
                dispersion[which(dispersion %in% colnames(mat))])
    }
  } else { 
    # used mcmc or vb
    if (x$stanfit@mode != 0) 
      STOP_no_draws()
    posterior <- rstan::extract(x$stanfit, permuted = FALSE, inc_warmup = FALSE)
    mat <- apply(posterior, 3L, FUN = function(y) y)
    if (!user_pars)
      pars <- grep("mean_PPD|log-posterior", # exclude these by default
                   colnames(mat), invert = TRUE, value = TRUE)
  }
  
  if (user_pars) {
    notfound <- which(!pars %in% colnames(mat))
    if (length(notfound)) 
      stop("No parameter(s) ", paste(pars[notfound], collapse = ", "), 
           call. = FALSE)
  }
  mat <- mat[, pars, drop = FALSE]
  if (!is.mer(x))
    return(mat)
  
  unpad_reTrms(mat)
}

#' @rdname as.matrix.stanreg
#' @export
as.data.frame.stanreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  mat <- as.matrix.stanreg(x, pars = pars, regex_pars = regex_pars, ...)
  as.data.frame(mat)
}
