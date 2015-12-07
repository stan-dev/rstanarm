#' Bayesian generalized linear additive models with group-specific terms via
#' Stan
#' 
#' Bayesian inference for GAMs with flexible priors.
#' 
#' @export
#' @templateVar fun stan_glmm4
#' @templateVar pkg gamm4
#' @templateVar pkgfun gamm4
#' @template return-stanreg-object
#' @template see-also
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' 
#' @param formula,random,family,data,knots,drop.unused.levels Same as for 
#'   \code{\link[gamm4]{gamm4}}.
#' @param subset,weights,na.action Same as \code{\link[stats]{glm}}, 
#'   but rarely specified.
#' @param ... Further arguments passed to \code{\link[rstan]{sampling}} 
#'   (e.g. \code{iter}, \code{chains},
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_gamm4} function is similar in syntax to 
#'   \code{\link[gamm4]{gamm4}}, which accepts a syntax that is similar to
#'   (but not quite as extensive as) that for \code{\link[mgcv]{gamm}} and 
#'   converts it internally into the syntax accepted by \code{\link[lme4]{glmer}}
#'   Rather than performing (restricted) maximum likelihood estimation the
#'   \code{stan_gamm4} function utilizes MCMC to perform Bayesian estimation. 
#'   The Bayesian model adds independent priors on the common regression
#'   coefficients (in the same way as \code{\link{stan_glm}}) and priors on 
#'   the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters, including the smooths. See
#'   \code{\link[gamm4]{gamm4}} for more information about the model
#'   specicification and \code{\link{priors}} for more information about the
#'   priors.
#' @importFrom lme4 getME
#' @examples
#' # see example(gamm4, package = "gamm4") but prefix gamm4() calls with stan_

stan_gamm4 <- function(formula, random = NULL, family = gaussian(), data = list(), 
                       weights = NULL, subset = NULL, na.action, knots = NULL, 
                       drop.unused.levels = TRUE, ..., 
                       prior = normal(), prior_intercept = normal(),
                       prior_ops = prior_options(),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "optimizing", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE) {

  mc <- match.call(expand.dots = FALSE)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame(2))
  if (is.function(family)) 
    family <- family()
  if (!is(family, "family"))
    stop("'family' must be a family")
  mc[[1]] <- quote(gamm4::gamm4)
  mc$... <- mc$prior <- mc$prior_intercept <- mc$prior_ops <- mc$prior_covariance <-
    mc$prior_PD <- mc$algorithm <- mc$adapt_delta <- NULL
  result <- suppressWarnings(eval(mc, parent.frame(1L)))
  group <- getME(result$mer, c("Zt", "cnms", "flist"))
  glmod <- getME(result$mer, c("X", "y"))               
  X <- glmod$X
  y <- glmod$y
  glmod$y <- NULL
  glmod$reTrms <- group

  if (TRUE) offset <- double(0)  
  else offset <- eval(attr(glmod$fr, "offset"), parent.frame(1L)) %ORifNULL% double(0)
  if (missing(weights)) weights <- double(0)
  if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
  if (is.null(prior)) prior <- list()
  if (is.null(prior_intercept)) prior_intercept <- list()
  if (!length(prior_ops)) prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)

  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, ...)
  call <- match.call(expand.dots = TRUE)
  prior.info <- get_prior_info(call, formals())
  
  Z <- pad_reTrms(Z = t(as.matrix(group$Zt)), cnms = group$cnms, 
                  flist = group$flist)$Z
  if (algorithm == "optimizing") 
    colnames(Z) <- .bnames(names(stanfit$par), value = TRUE)
  else colnames(Z) <- .bnames(names(stanfit), value = TRUE)
  
  fit <- nlist(stanfit, family, formula, offset, weights, x = cbind(X, Z), 
               y = y, data, prior.info, call, algorithm, glmod) 
  out <- stanreg(fit)
  # FIXME: replace guts of gam with point estimates from stanfit
  out$gam <- result$gam
  class(out) <- c(class(out), "lmerMod")
  return(out)
}
