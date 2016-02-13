# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015 Trustees of Columbia University
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

#' Bayesian generalized linear models via Stan
#'
#' Generalized linear modeling with optional prior distributions for 
#' the coefficients, intercept, and nuisance parameter.
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
#'   performed (if \code{algorithm} is \code{"sampling"}) via MCMC. The Bayesian
#'   model adds independent priors on the coefficients of the GLM. The 
#'   \code{stan_glm} function calls the workhorse \code{stan_glm.fit} function, 
#'   but it is also possible to call the latter directly.
#'   
#'   The \code{stan_glm.nb} function, which takes the extra argument
#'   \code{link}, is a simple wrapper for \code{stan_glm} with \code{family =
#'   \link{neg_binomial_2}(link)}.
#'   
#' @seealso The various vignettes for \code{stan_glm}.
#' 
#' @examples
#' if (!grepl("^sparc",  R.version$platform)) {
#' ### Linear regression
#' fit <- stan_glm(mpg / 10 ~ ., data = mtcars, QR = TRUE,
#'                 algorithm = "fullrank") # for speed only
#' plot(fit, ci_level = 0.5)
#' plot(fit, ci_level = 0.5, pars = "beta")
#' 
#' ### Logistic regression
#' data(lalonde, package = "arm")
#' dat <- within(lalonde, {
#'  re74_1k <- re74 / 1000
#'  re75_1k <- re75 / 1000
#' })
#' t7 <- student_t(df = 7)
#' fmla <- treat ~ re74_1k + re75_1k + educ + black + hisp + 
#'                married + nodegr + u74 + u75
#' fit2 <- stan_glm(fmla, data = dat, family = binomial(link="logit"), 
#'                  prior = t7, prior_intercept = t7, 
#'                  algorithm = "fullrank") # for speed only
#' plot(fit2, pars = c("black", "hisp", "nodegr", "u74", "u75"), 
#'      ci_level = 0.67, outer_level = 1, show_density = TRUE)
#' pp_check(fit2, check = "resid")
#' pp_check(fit2, check = "test", test = "mean")
#' }
#' \dontrun{
#' ### Poisson regression (example from help("glm")) 
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' fit3 <- stan_glm(counts ~ outcome + treatment, family = poisson(link="log"),
#'                  prior = normal(0, 1), prior_intercept = normal(0, 5))
#' plot(fit3, fill_color = "skyblue4", est_color = "maroon")
#' 
#' ### Gamma regression (example from help("glm"))
#' clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
#'                        lot1 = c(118,58,42,35,27,25,21,19,18),
#'                        lot2 = c(69,35,26,21,18,16,13,12,12))
#' fit4 <- stan_glm(lot1 ~ log_u, data = clotting, family = Gamma) 
#' print(fit4, digits = 2)                 
#' fit5 <- update(fit4, formula = lot2 ~ log_u)
#' }
#'
stan_glm <- function(formula, family = gaussian(), data, weights, subset,
                    na.action = NULL, offset = NULL, model = TRUE, 
                    x = FALSE, y = TRUE, contrasts = NULL, ..., 
                    prior = normal(), prior_intercept = normal(),
                    prior_ops = prior_options(), prior_PD = FALSE, 
                    algorithm = c("sampling", "optimizing", 
                                  "meanfield", "fullrank"),
                    adapt_delta = NULL, QR = FALSE) {
  
  algorithm <- match.arg(algorithm)
  family <- validate_family(family)
  if (missing(data)) 
    data <- environment(formula)
  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), 
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mf <- check_constant_vars(mf)
  mt <- attr(mf, "terms")
  Y <- array1D_check(model.response(mf, type = "any"))
  if (is.empty.model(mt))
    stop("No intercept or predictors specified.", call. = FALSE)
  X <- model.matrix(mt, mf, contrasts)
  weights <- validate_weights(as.vector(model.weights(mf)))
  offset <- validate_offset(as.vector(model.offset(mf)), y = Y)
  if (binom_y_prop(Y, family, weights)) {
    y1 <- as.integer(as.vector(Y) * weights)
    Y <- cbind(y1, y0 = weights - y1)
    weights <- double(0)
  }
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)

  stanfit <- stan_glm.fit(x = X, y = Y, weights = weights, 
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta, 
                          QR = QR, ...)
  fit <- nlist(stanfit, family, formula, offset, weights, x = X, y = Y, 
               data, prior.info = get_prior_info(call, formals()), 
               call = call, terms = mt, model = mf, 
               algorithm, na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"))
  out <- stanreg(fit)
  out$xlevels <- .getXlevels(mt, mf)
  if (!x) 
    out$x <- NULL
  if (!y) 
    out$y <- NULL
  if (!model) 
    out$model <- NULL
  
  return(out)
}

#' @rdname stan_glm
#' @export
#' @param link For \code{stan_glm.nb} only, the link function to use. See 
#'   \code{\link{neg_binomial_2}}.
stan_glm.nb <- function(..., link = "log") {
  if ("family" %in% names(list(...))) 
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call)) 
    names(call)[2L] <- "formula"
  mc[[1L]] <- quote(stan_glm)
  mc$link <- NULL
  mc$family <- neg_binomial_2(link = link)
  out <- eval(mc, parent.frame())
  out$call <- call
  
  return(out)
}
