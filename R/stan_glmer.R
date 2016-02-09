#' Bayesian generalized linear models with group-specific terms via Stan
#' 
#' Bayesian inference for GLMs with group-specific coefficients that have 
#' unknown covariance matrices with flexible priors.
#' 
#' @export
#' @templateVar armRef (Ch. 11-15)
#' @templateVar fun stan_glmer, stan_lmer, stan_glmer.nb
#' @templateVar pkg lme4
#' @templateVar pkgfun glmer
#' @template return-stanreg-object
#' @template see-also
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template reference-gelman-hill
#' 
#' @param formula,data,family Same as for \code{\link[lme4]{glmer}}.
#' @param subset,weights,offset Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param ... For \code{stan_glmer}, further arguments passed to 
#'   \code{\link[rstan]{sampling}} (e.g. \code{iter}, \code{chains}, 
#'   \code{cores}, etc.) or to \code{\link[rstan]{vb}} (if \code{algorithm} is 
#'   \code{"meanfield"} or \code{"fullrank"}). For \code{stan_lmer} and 
#'   \code{stan_glmer.nb}, \code{...} should also contain all relevant arguments
#'   to pass to \code{stan_glmer} (except \code{family}).
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_glmer} function is similar in syntax to 
#'   \code{\link[lme4]{glmer}} but rather than performing (restricted) maximum 
#'   likelihood estimation of generalized linear models, Bayesian estimation is 
#'   performed via MCMC. The Bayesian model adds independent priors on the 
#'   regression coefficients (in the same way as \code{\link{stan_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters. See \code{\link{priors}} for more information
#'   about the priors.
#'   
#'   The \code{stan_lmer} function is equivalent to \code{stan_glmer} with 
#'   \code{family = gaussian(link = "identity")}. 
#'   
#'   The \code{stan_glmer.nb} function, which takes the extra argument
#'   \code{link}, is a simple wrapper for \code{stan_glmer} with 
#'   \code{family = \link{neg_binomial_2}(link)}.
#'   
#'   
#' @seealso The vignette for \code{stan_glmer} and the \emph{Hierarchical 
#'   Partial Pooling} vignette.
#'    
#' @examples
#' # see help(example_model) for details on the model below
#' print(example_model, digits = 1)
#' 
#' @importFrom lme4 glFormula glmerControl
#' @importFrom Matrix Matrix t
stan_glmer <- function(formula, data = NULL, family = gaussian, 
                       subset, weights, 
                       na.action = getOption("na.action", "na.omit"), 
                       offset, contrasts = NULL, ...,
                       prior = normal(), prior_intercept = normal(),
                       prior_ops = prior_options(),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE) {
  
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  family <- validate_family(family)
  mc[[1]] <- quote(lme4::glFormula)
  mc$control <- glmerControl(check.nlev.gtreq.5 = "ignore",
                             check.nlev.gtr.1 = "stop",
                             check.nobs.vs.rankZ = "ignore",
                             check.nobs.vs.nlev = "ignore",
                             check.nobs.vs.nRE = "ignore")
  mc$prior <- mc$prior_intercept <- mc$prior_covariance <- mc$prior_ops <-
    mc$prior_PD <- mc$algorithm <- mc$scale <- mc$concentration <- mc$shape <-
    mc$adapt_delta <- mc$... <- mc$QR <- NULL
  glmod <- eval(mc, parent.frame())
  y <- glmod$fr[, as.character(glmod$formula[2L])]
  X <- glmod$X

  offset <- eval(attr(glmod$fr, "offset"), parent.frame()) %ORifNULL% double(0)
  weights <- validate_weights(weights)
  if (is.null(prior)) 
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, ...)

  Z <- pad_reTrms(Z = t(group$Zt), cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, x = cbind2(X, Z), 
               y = y, data, call, terms = NULL, model = NULL, 
               prior.info = get_prior_info(call, formals()),
               na.action, contrasts, algorithm, glmod)
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  
  return(out)
}

#' @rdname stan_glmer
#' @export
stan_lmer <- function(...) {
  if ("family" %in% names(list(...))) 
    stop("'family' should not be specified.")
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call)) 
    names(call)[2L] <- "formula"
  mc[[1L]] <- quote(stan_glmer)
  mc$REML <- NULL
  mc$family <- gaussian
  out <- eval(mc, parent.frame())
  out$call <- call
  
  return(out)
}


#' @rdname stan_glmer
#' @export
#' @param link For \code{stan_glmer.nb} only, the link function to use. See 
#'   \code{\link{neg_binomial_2}}.
#' 
stan_glmer.nb <- function(..., link = "log") {
  if ("family" %in% names(list(...))) 
    stop("'family' should not be specified.")
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call)) 
    names(call)[2L] <- "formula"
  mc[[1]] <- quote(stan_glmer)
  mc$REML <- mc$link <- NULL
  mc$family <- neg_binomial_2(link = link)
  out <- eval(mc, parent.frame())
  out$call <- call
  
  return(out)
}
