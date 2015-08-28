#' Fitting Bayesian generalized linear models with group-specific terms via
#' Stan
#'
#' Full Bayesian inference or optimization for generalized linear modeling with
#' group-specific terms with Gaussian, Student t, or Cauchy prior distributions
#' for the coefficients and flexible priors for the unknown covariance matrices.
#' 
#' @export
#' 
#' @templateVar fun stan_glmer, stan_lmer
#' @templateVar pkg lme4
#' @templateVar pkgfun glmer
#' @template return-stanreg-object
#' @template see-also
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-dots
#' 
#' @param formula,data,family Same as for \code{\link[lme4]{glmer}}.
#' @param subset,weights,offset Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely
#'   specified.
#' @param prior.for.covariance Cannot be \code{NULL}; see 
#'   \code{\link{decov}} for more information about the default arguments.   
#'
#' @details The \code{stan_glmer} function is similar in syntax to 
#'   \code{\link[lme4]{glmer}} but rather than performing (restricted) maximum 
#'   likelihood estimation of generalized linear models, full Bayesian 
#'   estimation is performed via Markov Chain Monte Carlo. The Bayesian model 
#'   adds independent Gaussian, Student t, or Cauchy priors on the coefficients 
#'   of the generalized linear model and priors on the terms of a decomposion
#'   of the covariance matrices of the group-specific parameters. See
#'   \code{\link{priors}} for more information about the priors.
#'   
#' @examples
#' # algorithm = "meanfield" is only for time constraints on examples
#' stan_glmer(mpg ~ . + (1|gear), data = mtcars, family = gaussian(),
#'            algorithm = "meanfield", tol_rel_obj = 0.05, seed = 12345)
#'
#' @importFrom lme4 glFormula
#' 
stan_glmer <- function (formula, data = NULL, family = gaussian, 
                        subset, weights, 
                        na.action = getOption("na.action", "na.omit"), 
                        offset, contrasts = NULL, ...,
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
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
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

  Z <- t(as.matrix(group$Zt))
  if (algorithm == "optimizing") 
    colnames(Z) <- grep("^b\\[", names(stanfit$par), value = TRUE)
  else colnames(Z) <- grep("^b\\[", names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, x = cbind(X, Z), 
               y = y, data, prior.info, call = match.call(expand.dots = TRUE), 
               terms = NULL, model = NULL, na.action, contrasts, algorithm)
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  
  out$glmod <- glmod # useful for post-processing
  
  return(out)
}

#' @rdname stan_glmer
#' @export
stan_lmer <- function (formula, data = NULL, subset, weights, na.action, offset, 
                       contrasts = NULL, ...,
                       prior = normal(), prior.for.intercept = normal(),
                       prior.options = prior_options(), 
                       prior.for.covariance = decov(), prior_PD = FALSE,
                      algorithm = c("sampling", "optimizing", "meanfield", 
                                    "fullrank")) {
  
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- quote(stan_glmer)
  mc$REML <- NULL
  mc$family <- gaussian
  return(eval(mc, parent.frame(1L)))
}
