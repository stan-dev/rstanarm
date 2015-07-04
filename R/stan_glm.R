# This file is part of rstanarm.
# Copyright 2013 Stan Development Team
# rstanarm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rstanarm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.


#' Fitting Bayesian generalized linear models via Stan
#'
#' Full Bayesian inference for generalized linear modeling with
#' Gaussian, Student t, or Cauchy prior distributions for the coefficients.
#'
#' @export
#'
#' @param formula,family,data,subset,x,y,model Same as 
#'   \code{\link[stats]{glm}}.
#' @param na.action,weights,offset,contrasts Same as 
#'   \code{\link[stats]{glm}}.
#' @param start Same as \code{\link[stats]{glm}}, but if not \code{NULL} also
#'   used as starting values for the MCMC. If \code{NULL} (the default), then
#'   \code{\link[rstan]{stan}} is initialized with \code{init = 'random'}.
#'   
#' @param prior Prior for coefficients. See \code{\link{priors}}.
#' @param prior.for.intercept Prior for intercept. See \code{\link{priors}}.
#' @param prior.options Additional options related to prior distributions. See
#'   \code{\link{priors}}.
#'   
#' @param ... Further arguments passed to \code{\link[rstan]{stan}} (e.g.
#'   \code{iter}, \code{chains}, \code{refresh}, etc.)
#' 
#'
#' @details The \code{stan_glm} function is similar in syntax to 
#'   \code{\link[stats]{glm}} but rather than performing maximum likelihood 
#'   estimation of generalized linear models, full Bayesian estimation is 
#'   performed via Markov Chain Monte Carlo. The Bayesian model adds independent
#'   Gaussian, Student t, or Cauchy priors on the coefficients of the 
#'   generalized linear model. The \code{stan_lm} function calls \code{stan_glm}
#'   with \code{family = gaussian}.
#' 
#' @return A list containing the components
#' \describe{
#'   \item{coefficients}{Named vector of coefficients (posterior means)}
#'   \item{residuals}{The residuals. For linear models \code{residuals}
#'    contains the response minus fitted values. Otherwise \code{residuals}
#'    contains the deviance residuals. See also \code{\link{residuals.stanreg}}}.
#'   \item{fitted.values}{The fitted mean values (for GLMs 
#'   the linear predictors are transformed by the invserse link function).}
#'   \item{linear.predictors}{The linear fit on the link scale (for linear models
#'   this is the same as \code{fitted.values}).}
#'   \item{covmat}{Variance-covariance matrix for the coefficients (estimated
#'   from the posterior draws.)}
#'   \item{y}{If requested, the \code{y} vector used.}
#'   \item{x}{If requested, the model matrix.}
#'   \item{model}{If requested, the model frame.}
#'   \item{family}{The \code{\link[stats]{family}} object used.}
#'   \item{prior.weights}{The weights provided by the user or a vector of
#'   \code{1}s if none were supplied.}
#'   \item{weights}{Same as \code{prior.weights} (note: this differs from
#'   \code{\link[stats]{glm}}, for which \code{prior.weights} and \code{weights} 
#'   are \emph{not} the same)}
#'   \item{df.residual}{The residual degrees of freedom}
#'   \item{call}{The matched call}
#'   \item{formula}{The formula supplied}
#'   \item{data}{The \code{data} argument}
#'   \item{prior.info}{A list with information about the prior distributions
#'   used}
#'   \item{stanfit}{The stanfit object returned by \code{\link[rstan]{stan}}}
#' } 
#' 
#' The accessor functions \code{coef}, \code{fitted}, and \code{resid}
#' can be used with objects of class \code{"stanreg"}. There are also  
#' \code{vcov}, \code{confint} and \code{\link{se}} methods.
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{lm}}, 
#' \code{\link[rstan]{stan}}
#' 
#' @examples 
#' \dontrun{
#' stan_lm(mpg ~ wt, data = mtcars)
#' stan_glm(mpg ~ wt, data = mtcars)
#' 
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' fit <- stan_lm(weight ~ group)
#' coef(fit)
#' 
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' stan_glm(counts ~ outcome + treatment, family = poisson())
#' }
#'

stan_glm <- function(formula, family = gaussian(), data, weights, subset,
                    na.action = NULL, start = NULL, offset = NULL, 
                    model = TRUE, x = FALSE, y = TRUE, contrasts = NULL,
                    prior = normal(), prior.for.intercept = normal(),
                    prior.options = prior_options(), 
                    ...) { # further arguments to stan()

  # Parse like glm()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if(!is(family, "family")) stop("'family' must be a family")
  
  if (missing(data)) data <- environment(formula)
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  arg_nms <- c("formula", "data", "subset", "weights", "na.action", "offset")
  m <- match(arg_nms, names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) names(Y) <- nm
  }

  if (!is.empty.model(mt)) X <- model.matrix(mt, mf, contrasts)
  else X <- matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
  if (is.null(weights)) weights <- rep(1.0, NROW(Y))
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  else offset <- rep(0, nrow(X))
  
  stanfit <- stan_glm.fit(x = X, y = Y, weights = weights, start = start, 
                          offset = offset, family = family, 
                          prior.dist = prior$dist,
                          prior.dist.for.intercept = prior.for.intercept$dist,
                          prior.mean = prior$location, prior.scale = prior$scale, 
                          prior.df = na_replace(prior$df, 1), 
                          prior.mean.for.intercept = prior.for.intercept$location, 
                          prior.scale.for.intercept = prior.for.intercept$scale,
                          prior.df.for.intercept = na_replace(prior.for.intercept$df, 1), 
                          scaled = prior.options$scaled, 
                          min.prior.scale = prior.options$min.prior.scale, 
                          prior.scale.for.dispersion = prior.options$prior.scale.for.dispersion, 
                          ...)
  
  # list of all the arguments and their values including any defaults (match.call
  # doesn't include defaults)
  all_args <- mget(names(formals()), sys.frame(sys.nframe()))
  prior.info <- all_args[grep("prior", names(all_args), fixed = TRUE)]
  
  fit <- nlist(stanfit, family, formula, offset, weights, x = X, y = Y, 
               data, prior.info, call = call, terms = mt, 
               model = if (model) mf else NULL, 
               na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"))
  
  fit <- stanreg(fit)
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL
  fit
}


