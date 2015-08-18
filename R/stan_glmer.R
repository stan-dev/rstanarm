#' Fitting Bayesian generalized linear models with group-specific terms via
#' Stan
#'
#' Full Bayesian inference or optimization for generalized linear modeling with
#' group-specific terms with Gaussian, Student t, or Cauchy prior distributions
#' for the coefficients and flexible priors for the unknown covariance matrices.
#' 
#' @export
#' 
#' @template return-stanreg-object
#' @templateVar fun stan_glmer, stan_lmer
#' @template see-also
#' @templateVar pkg MASS
#' @templateVar pkgfun polr
#' 
#' @param formula,data,family Same as for \code{\link[lme4]{glmer}}.
#' @param control,verbose,nAGQ,mustart,etastart,devFunOnly Same as for 
#'   \code{\link[lme4]{glmer}} but ignored.
#' @param subset,weights,na.action,offset,contrasts Same as 
#'   \code{\link[stats]{glm}}.
#' @param start If \code{NULL} (the default), then
#'   \code{\link[rstan]{stan}} is initialized with \code{init = 'random'}.
#'   If not \code{NULL} also used as starting values for the MCMC.
#' @param ... Further arguments passed to the function in the \pkg{rstan} 
#'   package named by \code{algorithm} (e.g., for the case of
#'   \code{\link[rstan]{sampling}}, \code{iter}, \code{chains}, etc.).
#' @param prior Prior for coefficients. Can be \code{NULL} to omit a prior
#'   and see \code{\link{priors}} otherwise.
#' @param prior.for.intercept Prior for intercept. Can be \code{NULL} to omit
#'   a prior and see \code{\link{priors}} otherwise.
#'   Note: The prior distribution for the intercept is set so it applies to
#'   the value when all common predictors are centered.
#' @param prior.options Additional options related to prior distributions. 
#'   Can be \code{NULL} to omit a prior on the dispersion and see
#'   \code{\link{priors}} otherwise.
#' @param prior.for.covariance Cannot be \code{NULL}; see 
#'   \code{\link{decov}} for more information about the default arguments.   
#' @param prior_PD A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to draw from the prior predictive distribution instead of
#'   conditioning on the outcome.
#' @param algorithm Character string (possibly abbreviated) among 
#'   \code{"sampling"}, \code{"optimizing"}, \code{"meanfield"}, and 
#'   \code{"fullrank"} indicating the estimation approach to use.
#'
#' @details The \code{stan_glmer} function is similar in syntax to 
#'   \code{\link[lme4]{glmer}} but rather than performing (restricted) maximum 
#'   likelihood estimation of generalized linear models, full Bayesian 
#'   estimation is performed via Markov Chain Monte Carlo. The Bayesian model 
#'   adds independent Gaussian, Student t, or Cauchy priors on the coefficients 
#'   of the generalized linear model and priors on the terms of a decomposion
#'   of the covariance matrices of the group-specific parameters. See
#'   \code{\link{priors}} for more information about the priors.
#' @examples
#' # algorithm = "meanfield" is only for time constraints on examples
#' stan_glmer(mpg ~ . + (1|gear), data = mtcars, family = gaussian(),
#'            algorithm = "meanfield", tol_rel_obj = 0.05, seed = 12345)
#'
#' @importFrom lme4 glFormula
#' 
stan_glmer <- function (formula, data = NULL, family = gaussian, 
                        control = NULL, start = NULL, verbose = 0L, 
                        nAGQ = 1L, subset, weights, 
                        na.action = getOption("na.action", "na.omit"), 
                        offset, contrasts = NULL, mustart, etastart, 
                        devFunOnly = FALSE, ...,
                        prior = normal(), prior.for.intercept = normal(),
                        prior.options = prior_options(),
                        prior.for.covariance = decov(), prior_PD = FALSE, 
                        algorithm = c("sampling", "optimizing", 
                                      "meanfield", "fullrank")) {
  
  mc <- match.call(expand.dots = FALSE)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame(2))
  if (is.function(family)) 
    family <- family()
  if(!is(family, "family"))
    stop("'family' must be a family")
  mc[[1]] <- quote(lme4::glFormula)
  mc$prior <- mc$prior.for.intercept <- mc$prior.options <- mc$prior_PD <-
    mc$algorithm <- mc$scale <- mc$concentration <- mc$shape <- mc$... <- NULL
  glmod <- eval(mc, parent.frame(1L))
  
  y <- glmod$fr[,as.character(glmod$formula[2])]
  X <- glmod$X

  if (missing(weights)) weights <- double(0)
  if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
  offset <- double(0)

  if (is.null(prior)) prior <- list()
  if (is.null(prior.for.intercept)) prior.for.intercept <- list()
  if (length(prior.options) == 0) {
    prior.options <- list(scaled = FALSE, prior.scale.for.dispersion = Inf)
  }
  group <- glmod$reTrms
  group$decov <- prior.for.covariance
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights, start = start, 
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
                          prior_PD = prior_PD, algorithm = algorithm, group = group, ...)
  
  all_args <- mget(names(formals()), sys.frame(sys.nframe()))
  prior.info <- all_args[grep("prior", names(all_args), fixed = TRUE)]
  mcout <- match.call(expand.dots = TRUE)

  fit <- nlist(stanfit, family, formula, offset, weights, x = cbind(X, group$Z), 
               y = y, data, prior.info, call = match.call(expand.dots = TRUE), 
               terms = NULL, model = NULL, na.action, contrasts, algorithm)
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  
  out$cnms <- glmod$reTrms$cnms # useful for post-processing
  
  return(out)
}

#' @rdname stan_glmer
#' @export
#' @param REML Ignored.
stan_lmer <- function (formula, data = NULL, REML = FALSE, control = NULL, 
                       start = NULL, verbose = 0L, subset, weights, na.action, offset, 
                       contrasts = NULL, devFunOnly = FALSE, 
                       prior = normal(), prior.for.intercept = normal(),
                       prior.options = prior_options(), 
                       prior.for.covariance = decov(), prior_PD = FALSE,
                      algorithm = c("sampling", "optimizing"),  ...) {
  
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- quote(stan_glmer)
  mc$REML <- NULL
  mc$family <- gaussian
  return(eval(mc, parent.frame(1L)))
}
