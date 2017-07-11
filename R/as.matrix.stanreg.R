# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

#' Extract the posterior sample
#' 
#' For models fit using MCMC (\code{algorithm="sampling"}), the posterior sample
#' ---the post-warmup draws from the posterior distribution--- can be extracted 
#' from a fitted model object as a matrix, data frame, or array. The 
#' \code{as.matrix} and \code{as.data.frame} methods merge all chains together, 
#' whereas the \code{as.array} method keeps the chains separate. For models fit 
#' using optimization (\code{"optimizing"}) or variational inference 
#' (\code{"meanfield"} or \code{"fullrank"}), there is no posterior sample but 
#' rather a matrix (or data frame) of 1000 draws from either the asymptotic
#' multivariate Gaussian sampling distribution of the parameters or the
#' variational approximation to the posterior distribution.
#' 
#' @method as.matrix stanreg
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template args-pars
#' @template args-regex-pars
#' @param ... Ignored.
#'   
#' @return A matrix, data.frame, or array, the dimensions of which depend on
#'   \code{pars} and \code{regex_pars}, as well as the model and estimation
#'   algorithm (see the Description section above).
#' 
#' @seealso \code{\link{stanreg-methods}}
#' 
#' @examples
#' \donttest{
#' if (!exists("example_model")) example(example_model)
#' # Extract posterior sample after MCMC
#' draws <- as.matrix(example_model)
#' print(dim(draws))
#' 
#' # For example, we can see that the median of the draws for the intercept 
#' # is the same as the point estimate rstanarm uses
#' print(median(draws[, "(Intercept)"]))
#' print(example_model$coefficients[["(Intercept)"]])
#' 
#' # The as.array method keeps the chains separate
#' draws_array <- as.array(example_model)
#' print(dim(draws_array)) # iterations x chains x parameters
#' 
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
  pars <- collect_pars(x, pars, regex_pars)
  user_pars <- !is.null(pars)
  
  if (used.optimizing(x)) {
    mat <- x$asymptotic_sampling_dist
    if (is.null(mat)) 
      STOP_no_draws()
    if (!user_pars) {
      aux <- c("sigma", "scale", "shape", "lambda", "reciprocal_dispersion")
      pars <- c(names(coef(x)), # return with coefficients first
                aux[which(aux %in% colnames(mat))])
    }
  } else { # used mcmc or vb
    mat <- as.matrix(x$stanfit)
    if (!user_pars)
      pars <- exclude_lp_and_ppd(colnames(mat))
  }
  if (user_pars)
    check_missing_pars(mat, pars)

  mat <- mat[, pars, drop = FALSE]
  if (!is.mer(x))
    return(mat)
  unpad_reTrms(mat)
}

#' @rdname as.matrix.stanreg
#' @method as.array stanreg
#' @export
as.array.stanreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  pars <- collect_pars(x, pars, regex_pars)
  if (!used.sampling(x))
    stop(
      "For models not fit using MCMC ", 
      "use 'as.matrix' instead of 'as.array'"
    )

  arr <- as.array(x$stanfit)
  if (identical(arr, numeric(0)))
    STOP_no_draws()
  
  if (!is.null(pars)) {
    check_missing_pars(arr, pars)
  } else {
    pars <- exclude_lp_and_ppd(last_dimnames(arr))
  }
  arr <- arr[, , pars, drop = FALSE]
  
  if (!is.mer(x))
    return(arr)
  unpad_reTrms(arr)
}


#' @rdname as.matrix.stanreg
#' @method as.data.frame stanreg
#' @export
as.data.frame.stanreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  mat <- as.matrix.stanreg(x, pars = pars, regex_pars = regex_pars, ...)
  as.data.frame(mat)
}



# internal ----------------------------------------------------------------
STOP_no_draws <- function() stop("No draws found.", call. = FALSE)

check_missing_pars <- function(x, pars) {
  notfound <- which(!pars %in% last_dimnames(x))
  if (length(notfound)) 
    stop(
      "No parameter(s) ", 
      paste(pars[notfound], collapse = ", "), 
      call. = FALSE
    )
}

exclude_lp_and_ppd <- function(pars) {
  grep(
    pattern = "mean_PPD|log-posterior", 
    x = pars, 
    invert = TRUE, 
    value = TRUE
  )
}

