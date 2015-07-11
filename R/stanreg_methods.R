#' Methods
#' @name stanreg_methods
#' 
#' @export
#' 
#' @param object,x a model fit with \code{\link{stan_lm}} or
#'   \code{\link{stan_glm}}.
#' @param ... other arguments. See Details.
#' @param parm a character vector of parameter names.
#' @param level confidence level.
#' @details For \code{residuals}, currently the only argument that can be 
#'   specified in \code{...} is \code{type}, which is used to indicate the type 
#'   of residuals to be  returned. \code{type} can be one of \code{'response'},
#'   \code{'deviance'}, or \code{'pearson'}. If omitted, the default is
#'   \code{'response'} for linear models and otherwise the default is
#'   \code{'deviance'}.
#' 
residuals.stanreg <- function(object, ...) {
  type <- c(...)
  nt <- is.null(type)
  if (!nt) {
    if (!(type %in% c("response", "deviance", "pearson")))
      stop("'type' should be 'response', 'deviance', or 'pearson'")
  }
  fam <- object$family
  if (fam$family == "gaussian" && fam$link == "identity") {
    if (nt) type <- "response"
    rr <- residuals.lm(object, type = type)
  } else {
    if (nt) type <- "deviance"
    rr <- residuals.glm(object, type = type)  
  } 
  attr(rr, "type") <- unname(type)
  rr
}

#' @rdname stanreg_methods
#' @export 
vcov.stanreg <- function(object, ...) {
  object$covmat
}

#' @rdname stanreg_methods
#' @export
se <- function(object, parm) UseMethod("se")

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
  print(x$stanfit, pars = "lp__", include = FALSE)
}

#' @rdname stanreg_methods
#' @export
summary.stanreg <- function(object, ...) {
  # use RStan's summary just as placeholder. we should replace this with our own
  # summary 
  summary(object$stanfit, ...)$summary
}

