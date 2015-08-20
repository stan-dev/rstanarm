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
#' Full Bayesian inference or optimization for generalized linear modeling with 
#' Gaussian, Student t, or Cauchy prior distributions for the coefficients.
#'
#' @export
#' 
#' @template return-stanreg-object
#' @templateVar fun stan_glm
#' @template return-stanfit-object
#' @templateVar fitfun stan_glm.fit
#' @template see-also
#' @templateVar pkg stats
#' @templateVar pkgfun glm
#' 
#'
#' @param formula,family,data,subset Same as \code{\link[stats]{glm}}.
#' @param x,y In \code{stan_glm}, logical scalars indicating whether to
#'   return the design matrix and response vector. In \code{stan_glm.fit},
#'   a design matrix and response vector.   
#' @param model,na.action,weights,offset,contrasts Same as 
#'   \code{\link[stats]{glm}}.
#' @param ... Further arguments passed to the function in the \pkg{rstan} 
#'   package named by \code{algorithm} (e.g., for the case of
#'   \code{\link[rstan]{sampling}}, \code{iter}, \code{chains}, etc.).
#' @param prior Prior for coefficients. Can be \code{NULL} to omit a prior
#'   and see \code{\link{priors}} otherwise.
#' @param prior.for.intercept Prior for intercept. Can be \code{NULL} to omit
#'   a prior and see \code{\link{priors}} otherwise.
#'   Note: The prior distribution for the intercept is set so it applies to
#'   the value when all predictors are centered.
#' @param prior.options Additional options related to prior distributions. 
#'   Can be \code{NULL} to omit a prior on the dispersion and see
#'   \code{\link{priors}} otherwise.
#' @param prior_PD A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to draw from the prior predictive distribution instead of
#'   conditioning on the outcome.
#' @param algorithm Character string (possibly abbreviated) among 
#'   \code{"sampling"}, \code{"optimizing"}, \code{"meanfield"}, and 
#'   \code{"fullrank"} indicating the estimation approach to use.
#' 
#' @details The \code{stan_glm} function is similar in syntax to 
#'   \code{\link[stats]{glm}} but rather than performing maximum likelihood 
#'   estimation of generalized linear models, full Bayesian estimation is 
#'   performed via Markov Chain Monte Carlo. The Bayesian model adds independent
#'   Gaussian, Student t, or Cauchy priors on the coefficients of the 
#'   generalized linear model. The \code{stan_glm} function calls the workhorse
#'   \code{stan_glm.fit} function, but it is possible to call the latter
#'   directly.
#' 
#' @examples 
#' # algorithm = "meanfield" is only for time constraints on examples
#' stan_glm(mpg ~ ., data = mtcars, algorithm = "meanfield", seed = 12345)
#'  
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' stan_glm(counts ~ outcome + treatment, family = poisson(), 
#'          algorithm = "meanfield")
#'

stan_glm <- function(formula, family = gaussian(), data, weights, subset,
                    na.action = NULL, offset = NULL, model = TRUE, 
                    x = FALSE, y = TRUE, contrasts = NULL, ..., 
                    prior = normal(), prior.for.intercept = normal(),
                    prior.options = prior_options(), 
                    prior_PD = FALSE,  algorithm = c("sampling", "optimizing", 
                                                     "meanfield", "fullrank")){

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
  else X <- matrix(NA_real_, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
  if (is.null(weights)) weights <- double(0) #rep(1.0, NROW(Y))
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  else offset <- double(0) #rep(0, nrow(X))
  
  # if Y is proportion of successes and weights is total number of trials
  if (family$family == "binomial" && NCOL(Y) == 1L && all(Y > 0 & Y <= 1)) {
      if (!identical(weights, double(0)) && all(weights > 0)) {
        y1 <- as.integer(as.vector(Y) * weights)
        Y <- cbind(y1, weights - y1)
        weights <- double(0)
    }
  }
  
  if (is.null(prior)) prior <- list()
  if (is.null(prior.for.intercept)) prior.for.intercept <- list()
  if (length(prior.options) == 0) {
    prior.options <- list(scaled = FALSE, prior.scale.for.dispersion = Inf)
  }
  algorithm <- match.arg(algorithm)

  stanfit <- stan_glm.fit(x = X, y = Y, weights = weights, 
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
                          prior_PD = prior_PD, algorithm = algorithm, ...)
  
  # list of all the arguments and their values including any defaults (match.call
  # doesn't include defaults)
  all_args <- mget(names(formals()), sys.frame(sys.nframe()))
  prior.info <- all_args[grep("prior", names(all_args), fixed = TRUE)]
  
  fit <- nlist(stanfit, family, formula, offset, weights, x = X, y = Y, 
               data, prior.info, call = call, terms = mt, model = mf, 
               algorithm, na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"))
  
  out <- stanreg(fit)
  if (!x) out$x <- NULL
  if (!y) out$y <- NULL
  if (!model) out$model <- NULL
  out
}


