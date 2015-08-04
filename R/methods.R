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
  if(x$stanfit@mode == 0) print(x$stanfit, pars = "lp__", include = FALSE, ...)
  else if (is.null(x$family)) {
    mark <- c(names(x$coefficients), 
              grep("|", rownames(x$stan_summary), fixed = TRUE, value = TRUE))
    print(x$stan_summary[mark,,drop = FALSE], ...)
  }
  else {
    mark <- names(x$coefficients)
    if (x$family$family == "gaussian") mark <- c(mark, "sigma")
    else if (x$family$family == "Negative Binomial") mark <- c(mark, "overdispersion")
    print(x$stan_summary[mark,,drop=FALSE], ...)
  }
}

#' @rdname stanreg_methods
#' @export
summary.stanreg <- function(object, ...) {
  # use RStan's summary just as placeholder. we should replace this with our own
  # summary 
  if(object$stanfit@mode == 0) summary(object$stanfit, ...)$summary
  else {
    mark <- names(object$coefficients)
    if (object$family$family == "gaussian") mark <- c(mark, "sigma")
    else if (object$family$family == "Negative Binomial") mark <- c(mark, "overdispersion")
    object$stan_summary[mark,,drop=FALSE]
  }
}

#' @rdname stanreg_methods
#' @export
log_lik.stanreg <- function(object) {
  object$log_lik
}

