#' Methods
#' @name stanreg_methods
#' 
#' @export
#' 
#' @param object,x a model fit with \code{\link{stan_lm}} or
#'   \code{\link{stan_glm}}.
#' @param ... other arguments to \code{print} or \code{summary}. See Details.
#' @param parm a character vector of parameter names.
#' @param level confidence level.
#' @note Unlike \code{\link[stats]{glm}}, residuals are of type \code{'response'} 
#' not \code{'deviance'} (see \code{\link[stats]{residuals.glm}}). 
#' 
residuals.stanreg <- function(object, ...) {
  object$residuals
}

#' @rdname stanreg_methods
#' @export 
vcov.stanreg <- function(object, ...) {
  object$covmat
}

#' @rdname stanreg_methods
#' @export
se.stanreg <- function(object, parm) {
  pnms <- names(coef(object))
  if (missing(parm)) parm <- pnms
  else if (is.numeric(parm)) parm <- pnms[parm]
  object$stan_summary[parm, "sd"]
}

#' @rdname stanreg_methods
#' @export
confint.stanreg <- function (object, parm, level = 0.95, ...) {
  # just a placeholder. we should replace this with a confint method that
  # returns posterior quantiles probably
  confint.default(object, parm, level, ...)
}

#' @rdname stanreg_methods
#' @export
print.stanreg <- function(x, ...) {
  # use RStan's print just as placeholder. we should replace this with our own
  # print method
  print(x$stanfit)
}

#' @rdname stanreg_methods
#' @export
summary.stanreg <- function(object, ...) {
  # use RStan's summary just as placeholder. we should replace this with our own
  # summary 
  summary(object$stanfit, ...)$summary
}

#' @rdname stanreg_methods
#' @export
log_lik.stanreg <- function(object) {
  object$log_lik
}

