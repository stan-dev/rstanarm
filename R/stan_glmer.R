#' Bayesian generalized linear models with group-specific terms via Stan
#'
#' Bayesian inference for GLMs with group-specific coefficients that have
#' unknown covariance matrices with flexible priors.
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
#' @param prior_covariance Cannot be \code{NULL}; see 
#'   \code{\link{decov}} for more information about the default arguments.   
#'
#' @details The \code{stan_glmer} function is similar in syntax to 
#'   \code{\link[lme4]{glmer}} but rather than performing (restricted) maximum
#'   likelihood estimation of generalized linear models, Bayesian estimation
#'   is performed (if \code{algorithm = "sampling"}) via MCMC. The Bayesian 
#'   model adds independent priors on the coefficients (in the same way as
#'   \code{\link{stan_glm}}) and priors on the terms of a decomposition
#'   of the covariance matrices of the group-specific parameters. See
#'   \code{\link{priors}} for more information about the priors.
#'   
#' @examples
#' \dontrun{ 
#' options(mc.cores = parallel::detectCores())
#' data(cake, package = "lme4")
#' stan_glmer(angle ~ recipe + temp + (1|recipe:replicate), data = cake, seed = 12345)
#' }
#' @importFrom lme4 glFormula
#' 
stan_glmer <- function (formula, data = NULL, family = gaussian, 
                        subset, weights, 
                        na.action = getOption("na.action", "na.omit"), 
                        offset, contrasts = NULL, ...,
                        prior = normal(), prior_intercept = normal(),
                        prior_ops = prior_options(),
                        prior_covariance = decov(), prior_PD = FALSE, 
                        algorithm = c("sampling", "optimizing")) {
  
  mc <- match.call(expand.dots = FALSE)
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame(2))
  if (is.function(family)) 
    family <- family()
  if (!is(family, "family"))
    stop("'family' must be a family")
  mc[[1]] <- quote(lme4::glFormula)
  mc$prior <- mc$prior_intercept <- mc$prior_ops <- mc$prior_PD <-
    mc$algorithm <- mc$scale <- mc$concentration <- mc$shape <- mc$... <- NULL
  glmod <- eval(mc, parent.frame(1L))
  y <- glmod$fr[,as.character(glmod$formula[2])]
  X <- glmod$X

  offset <- eval(attr(glmod$fr, "offset"), parent.frame(1L)) %ORifNULL% double(0)
  if (missing(weights)) weights <- double(0)
  if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
  if (is.null(prior)) prior <- list()
  if (is.null(prior_intercept)) prior_intercept <- list()
  if (!length(prior_ops)) prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, group = group, ...)
  all_args <- mget(names(formals()), sys.frame(sys.nframe()))
  prior.info <- all_args[grep("prior", names(all_args), fixed = TRUE)]

  Z <- t(as.matrix(group$Zt))
  if (algorithm == "optimizing") 
    colnames(Z) <- grep("^b\\[", names(stanfit$par), value = TRUE)
  else colnames(Z) <- grep("^b\\[", names(stanfit), value = TRUE)
  
  fit <- nlist(stanfit, family, formula, offset, weights, x = cbind(X, Z), 
               y = y, data, prior.info, call = match.call(expand.dots = TRUE), 
               terms = NULL, model = NULL, na.action, contrasts, algorithm)
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  out$glmod <- glmod
  return(out)
}

#' @rdname stan_glmer
#' @export
stan_lmer <- function (formula, data = NULL, subset, weights, na.action, offset, 
                       contrasts = NULL, ...,
                       prior = normal(), prior_intercept = normal(),
                       prior_ops = prior_options(), 
                       prior_covariance = decov(), prior_PD = FALSE,
                       algorithm = c("sampling", "optimizing")) {
  
  mc <- call <- match.call(expand.dots = TRUE)
  mc[[1]] <- quote(stan_glmer)
  mc$REML <- NULL
  mc$family <- gaussian
  out <- eval(mc, parent.frame(1L))
  out$call <- call
  return(out)
}
