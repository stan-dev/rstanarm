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


#' Bayesian generalized linear models via Stan
#'
#' Generalized linear modeling with optional prior distributions for 
#' the coefficients, intercept, and nuisance parameter
#'
#' @export
#' @templateVar armRef (Ch. 3-6)
#' @templateVar pkg stats
#' @templateVar pkgfun glm
#' @templateVar sameargs model,offset,weights 
#' @templateVar rareargs na.action,contrasts
#' @templateVar fun stan_glm, stan_glm.nb
#' @templateVar fitfun stan_glm.fit
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template see-also
#' @template args-formula-data-subset
#' @template args-same-as
#' @template args-same-as-rarely
#' @template args-x-y
#' @template args-dots
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template reference-gelman-hill
#' 
#' @param family Same as \code{\link[stats]{glm}}, except negative binomial GLMs
#'   are also possible using the \code{\link{neg_binomial_2}} family object.
#' 
#' @details The \code{stan_glm} function is similar in syntax to 
#'   \code{\link[stats]{glm}} but rather than performing maximum likelihood 
#'   estimation of generalized linear models, full Bayesian estimation is 
#'   performed (if \code{algorithm = "sampling"}) via MCMC. The Bayesian model 
#'   adds independent priors on the coefficients of the GLM. The 
#'   \code{stan_glm} function calls the workhorse \code{stan_glm.fit} function,
#'   but it is possible to call the latter directly.
#'   
#'   The \code{stan_glm.nb} function, which takes the extra argument
#'   \code{link}, is a simple wrapper for \code{stan_glm} with \code{family =
#'   \link{neg_binomial_2}(link)}.
#' 
#' @examples 
#' \dontrun{ 
#' ### Linear regression
#' fit <- stan_glm(mpg ~ ., data = mtcars, seed = 12345, cores = 1)
#' plot(fit, ci_level = 0.5)
#' 
#' ### Poisson regression  
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' fit2 <- stan_glm(counts ~ outcome + treatment, family = poisson(link="log"),
#'                  prior = normal(0, 2.5), prior_intercept = normal(0, 10), 
#'                  cores = 1)
#' plot(fit2, ci_level = 0.75, outer_level = 0.99, show_density = TRUE)
#' 
#' ### Logistic regression
#' data(lalonde, package = "arm")
#' ?lalonde
#' t7 <- student_t(df = 7) 
#' f <- treat ~ re74 + re75 + educ + black + hisp + married + nodegr + u74 + u75
#' fit3 <- stan_glm(f, data = lalonde, family = binomial(link="logit"), 
#'                  prior = t7, prior_intercept = t7, cores = 1)
#' plot(fit3)
#' ppcheck(fit3, check = "resid")
#' ppcheck(fit3, check = "test", test = "mean")
#' }
#'
stan_glm <- function(formula, family = gaussian(), data, weights, subset,
                    na.action = NULL, offset = NULL, model = TRUE, 
                    x = FALSE, y = TRUE, contrasts = NULL, ..., 
                    prior = normal(), prior_intercept = normal(),
                    prior_ops = prior_options(), prior_PD = FALSE, 
                    algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                    adapt_delta = NULL, QR = FALSE) {

  # Parse like glm()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (!is(family, "family")) stop("'family' must be a family")
  
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
  
  # check if any x variables are constants 
  is_constant <- apply(mf, 2, FUN = function(x) length(unique(x))) == 0
  if (any(is_constant)) stop("Constant variable(s) found: ",
                             paste(colnames(mf)[is_constant], collapse = ", "))

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
  if (is.binomial(family$family) && NCOL(Y) == 1L && is.numeric(Y)) { 
    if (all(findInterval(Y, c(.Machine$double.eps,1)) == 1)) { 
      if (!identical(weights, double(0)) && all(weights > 0)) {
        y1 <- as.integer(as.vector(Y) * weights)
        Y <- cbind(y1, weights - y1)
        weights <- double(0)
      }
    }
  }
  
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  algorithm <- match.arg(algorithm)

  stanfit <- stan_glm.fit(x = X, y = Y, weights = weights, 
                          offset = offset, family = family,
                          prior = prior,
                          prior_intercept = prior_intercept,
                          prior_ops = prior_ops,
                          prior_PD = prior_PD, algorithm = algorithm, 
                          adapt_delta = adapt_delta, QR = QR, ...)
  fit <- nlist(stanfit, family, formula, offset, weights, x = X, y = Y, 
               data, prior.info = get_prior_info(call, formals()), 
               call = call, terms = mt, model = mf, 
               algorithm, na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"))
  
  out <- stanreg(fit)
  if (!x) out$x <- NULL
  if (!y) out$y <- NULL
  if (!model) out$model <- NULL
  out
}

#' @rdname stan_glm
#' @export
#' @param link For \code{stan_glm.nb} only, the link function to use. See 
#'   \code{\link{neg_binomial_2}}.
stan_glm.nb <- function(..., link = "log") {
  mc <- call <- match.call()
  mc[[1L]] <- quote(stan_glm)
  mc$link <- NULL
  mc$family <- neg_binomial_2(link = link)
  out <- eval(mc, parent.frame(1L))
  out$call <- call
  out
}
