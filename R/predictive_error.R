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

#' In-sample or out-of-sample predictive errors
#' 
#' This is a convenience function for computing \eqn{y - y^{rep}}{y - yrep}
#' (in-sample) or \eqn{y - \tilde{y}}{y - ytilde} (out-of-sample) for each draw
#' from the posterior predictive distribution.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param newdata,draws,seed,offset,re.form Optional arguments passed to 
#'   \code{\link{posterior_predict}}. For binomial models, please see the
#'   \strong{Note} section below if \code{newdata} will be specified.
#' @param ... Currently ignored by the method for stanreg objects. The S3
#'   generic uses \code{...} to pass arguments to any defined methods.
#' 
#' @return A \code{draws} by \code{nrow(newdata)} matrix. If \code{newdata} is 
#'   not specified then it will be \code{draws} by \code{nobs(object)}. Each row
#'   of the matrix is a vector of predictive errors.
#'   
#' @note The \strong{Note} section in \code{\link{posterior_predict}} about 
#'   \code{newdata} for binomial models also applies for
#'   \code{predictive_error}, with one important difference. For
#'   \code{posterior_predict} if the left-hand side of the model formula is 
#'   \code{cbind(successes, failures)} then the particular values of 
#'   \code{successes} and \code{failures} in \code{newdata} don't matter, only 
#'   that they add to the desired number of trials. \strong{This is not the case
#'   for} \code{predictive_error}. For \code{predictive_error} the particular
#'   value of \code{successes} matters because it is used as \eqn{y} when
#'   computing the error.
#' 
#' @seealso \code{\link[=posterior_predict.stanreg]{posterior_predict}} to draw
#'   from the posterior predictive distribution without computing predictive
#'   errors.
#'   
#' @examples
#' if (!exists("example_model")) example(example_model)
#' err1 <- predictive_error(example_model, draws = 50)
#' hist(err1)
#' 
#' # Using newdata with a binomial model
#' formula(example_model)
#' nd <- data.frame(
#'  size = c(10, 20), 
#'  incidence = c(5, 10), 
#'  period = factor(c(1,2)), 
#'  herd = c(1, 15)
#' )
#' err2 <- predictive_error(example_model, newdata = nd, draws = 10, seed = 1234)
#' 
#' 
#' fit <- stan_glm(mpg ~ wt, data = mtcars, iter = 300)
#' hist(predictive_error(fit))
#' 
predictive_error <- function(object, ...) {
  UseMethod("predictive_error")
}

#' @rdname predictive_error
#' @export 
predictive_error.stanreg <- function(object, 
                                     newdata = NULL, 
                                     draws = NULL, 
                                     re.form = NULL, 
                                     seed = NULL, 
                                     offset = NULL, ...) {
  if (used.optimizing(object))
    STOP_not_optimizing("predictive_error")
  if (inherits(object, "polr"))
    stop("'predictive_error' is not currently available for stan_polr.")
  
  y <- if (is.null(newdata))
    get_y(object) else eval(formula(object)[[2L]], newdata)
  
  fam <- family(object)$family
  if (is.binomial(fam) && NCOL(y) == 2)
    y <- y[, 1]
  
  preds <- posterior_predict(
    object,
    newdata = newdata,
    draws = draws,
    offset = offset,
    seed = seed,
    re.form = re.form,
    ...
  )
  err <- sweep(preds, MARGIN = 2L, STATS = as.array(y), FUN = "-")
  as.matrix(-1 * err)
}
