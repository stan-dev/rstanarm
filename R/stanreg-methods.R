#' Methods
#' 
#' Methods for \link[=stanreg-objects]{'stanreg' objects}.
#' 
#' @name stanreg-methods
#' 
#' @param object,x A fitted model object returned by one of the \pkg{rstanarm} 
#'   modeling functions. This will be a list with class 'stanreg' as well as at 
#'   least one of 'lm', 'glm', 'polr', 'lmerMod', or 'aov'.
#' @param ... Ignored.
#' @param parm A character vector of parameter names.
#' @param level The confidence level to use.
#' @details The \code{se} method returns standard errors and the \code{log_lik} 
#'   method returns the pointwise log-likelihood matrix. Unlike 
#'   \code{\link[stats]{residuals.glm}}, residuals are of type \code{'response'}
#'   not \code{'deviance'}.
#'
#' @seealso \code{\link{stanreg-objects}}
#' 
NULL

#' @rdname stanreg-methods
#' @export 
residuals.stanreg <- function(object, ...) {
  object$residuals
}

#' @rdname stanreg-methods
#' @export 
vcov.stanreg <- function(object, ...) {
  object$covmat
}

#' @rdname stanreg-methods
#' @export
confint.stanreg <- function (object, parm, level = 0.95, ...) {
  # just a placeholder. we should replace this with a confint method that
  # returns posterior quantiles probably
  confint.default(object, parm, level, ...)
}

#' @rdname stanreg-methods
#' @export
coef.stanreg <- function(object, ...)  {
  object$coefficients
}

#' @rdname stanreg-methods
#' @export
fitted.stanreg <- function(object, ...)  {
  object$fitted.values
}

# standard errors
se <- function(object, parm) UseMethod("se")
#' @rdname stanreg-methods
#' @export
se.stanreg <- function(object, parm) {
  pnms <- names(coef(object))
  if (missing(parm)) parm <- pnms
  else if (is.numeric(parm)) parm <- pnms[parm]
  object$stan_summary[parm, "sd"]
}

# Compute pointwise log-likelihood matrix
log_lik <- function(object, ...) UseMethod("log_lik")
#' @rdname stanreg-methods
#' @export
log_lik.stanreg <- function(object, ...) {
  if (object$algorithm != "sampling")
    stop("Only available for MCMC.", call. = FALSE)
  fun <- .llfun(object)
  args <- .llargs(object)
  sapply(seq_len(args$N), function(i) {
    as.vector(fun(i = i, data = args$data, draws = args$draws)) 
  })
}
