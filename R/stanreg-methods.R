#' Methods
#' 
#' @name stanreg-methods
#' 
#' @param object,x A fitted model object returned by one of the \pkg{rstanarm} 
#'   modeling functions. This will be a list with class 'stanreg' as well as at 
#'   least one of 'lm', 'glm', 'polr', 'lmerMod', or 'aov'.
#' @param ... Other arguments to \code{print} or \code{summary}. See Details.
#' @param parm A character vector of parameter names.
#' @param level The confidence level to use.
#' @note Unlike \code{\link[stats]{glm}}, residuals are of type
#'   \code{'response'} not \code{'deviance'} (see
#'   \code{\link[stats]{residuals.glm}}).
#'
#' @seealso \code{\link{stanreg-objects}}
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
se.stanreg <- function(object, parm) {
  pnms <- names(coef(object))
  if (missing(parm)) parm <- pnms
  else if (is.numeric(parm)) parm <- pnms[parm]
  object$stan_summary[parm, "sd"]
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

#' @rdname stanreg-methods
#' @export
log_lik.stanreg <- function(object) {
  object$log_lik
}

#' @rdname stanreg-methods
#' @export
print.stanreg <- function(x, ...) {
  # use RStan's print just as placeholder. we should replace this with our own
  # print method
  if (x$stanfit@mode == 0) print(x$stanfit, pars = "lp__", include = FALSE, ...)
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
  if (is(x, "aov")) {
    labels <- attributes(x$terms)$term.labels
    patterns <- gsub(":", ".*:", labels)
    groups <- sapply(patterns, simplify = FALSE, FUN = grep, 
                     x = dimnames(extract(x$stanfit, pars = "beta", permuted = FALSE))$parameters)
    names(groups) <- gsub(".*", "", names(groups), fixed = TRUE)
    groups <- groups[sapply(groups, length) > 0]
    effects_dim <- dim(x$effects)
    effects <- x$effects^2
    effects <- sapply(groups, FUN = function(i) apply(effects[,,i,drop = FALSE], 1:2, mean))
    dim(effects) <- c(effects_dim[-3], ncol(effects))
    dimnames(effects) <- list(iterations = NULL, chains = paste0("chain", 1:ncol(x$stanfit), sep = ":"),
                              parameters = paste("Mean Sq", names(groups)))
    cat("ANOVA-like table\n")
    rstan::monitor(effects, warmup = 0, print = TRUE, digits_summary = 2)
  }
}

#' @rdname stanreg-methods
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
