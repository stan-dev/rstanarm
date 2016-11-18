#' Extract and/or print a summary of the priors used for an rstanarm model
#'
#' @aliases prior_summary
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param digits Number of digits to use for rounding.
#' @param ... Currently ignored by the method for stanreg objects. The S3
#'   generic uses \code{...} to pass arguments to any defined methods.
#' 
#' @details For some models you may see "\code{adjusted scale}" in the printed 
#'   output and adjusted scales included in the object returned by 
#'   \code{prior_summary}. These adjusted scale values are the prior scales 
#'   actually used by \pkg{rstanarm} and are computed by adjusting the prior
#'   scales specified by the user to account for the scales of the predictors
#'   (as described in the documentation for the \code{scaled} argument to 
#'   \code{\link{prior_options}}). For models with adjusted prior scales, 
#'   refitting the model with \code{prior_ops=prior_options(scaled=FALSE)} will
#'   disable this feature.
#' @return A list of class "prior_summary.stanreg", which has its own print
#'   method.
#'   
#' @examples
#' if (!exists("example_model")) example(example_model) 
#' prior_summary(example_model)
#' 
#' priors <- prior_summary(example_model)
#' names(priors)
#' priors$prior$scale
#' priors$prior$adjusted_scale
#' 
#' # for a glm with adjusted scales (see Details, above), compare 
#' # the default (rstanarm adjusting the scales) to setting 
#' # prior_ops=prior_options(scaled=FALSE)
#' fit <- stan_glm(mpg ~ wt + am, data = mtcars, 
#'                 prior = normal(0, c(2.5, 4)), 
#'                 prior_intercept = normal(0, 5), 
#'                 iter = 10, chains = 1) # only for demonstration 
#' prior_summary(fit)
#' 
#' fit2 <- update(fit, prior_ops = prior_options(scaled = FALSE))
#' prior_summary(fit2)
#' 
prior_summary.stanreg <- function(object, digits = 2,...) {
  x <- object[["prior.info"]]
  if (is.null(x)) {
    message("Priors not found in stanreg object.")
    return(NULL)
  }
  structure(x, class = "prior_summary.stanreg", 
            model_name = deparse(substitute(object)), 
            print_digits = digits)
}

#' @export
#' @method print prior_summary.stanreg
print.prior_summary.stanreg <- function(x, digits, ...) {
  if (missing(digits))
    digits <- attr(x, "print_digits") %ORifNULL% 2
  .dig <- digits
  .fr2 <- function(y, .digits = .dig, ...) format(y, digits = .digits, ...)
  .fr3 <- function(y, .nsmall = .dig) .fr2(y, nsmall = .nsmall)
  
  prior_intercept <- x[["prior_intercept"]]
  prior_coef <- x[["prior"]]
  prior_covariance <- x[["prior_covariance"]]
  prior_counts <- x[["prior_counts"]] # for stan_polr
  
  msg <- paste0("Priors for model '", attr(x, "model_name"), "'")
  cat(msg, "\n------")
  
  if (!is.null(prior_intercept)) {
    p <- prior_intercept
    cat("\nIntercept\n ~",
        if (is.na(p$dist)) {
          "flat"
        } else if (is.null(p$df)) {
          paste0(p$dist,"(location = ", .fr2(p$location), 
                 ", scale = ", .fr2(p$scale),")")
        } else {
          paste0(p$dist, "(df = ", .fr2(p$df), ", 
                 location = ", .fr2(p$location), 
                 ", scale = ", .fr2(p$scale), ")")
        }
      )
    if (!is.null(p$adjusted_scale))
      cat("\n     **adjusted scale =", .fr3(p$adjusted_scale))
  }
  
  if (!is.null(prior_coef)) {
    p <- prior_coef
    if (!(p$dist %in% c("R2", NA))) {
      if (p$dist %in% c("normal", "student_t", "cauchy")) {
        p$location <- .format_pars(p$location, .fr2)
        p$scale <- .format_pars(p$scale, .fr2)
        if (!is.null(p$df))
          p$df <- .format_pars(p$df, .fr2)
        if (!is.null(p$adjusted_scale))
          p$adjusted_scale <- .format_pars(p$adjusted_scale, .fr2)
      } else if (p$dist %in% c("hs_plus")) {
        p$df1 <- .format_pars(p$df, .fr2)
        p$df2 <- .format_pars(p$scale, .fr2)
      } else if (p$dist %in% c("hs")) {
        p$df <- .format_pars(p$df, .fr2)
      }
    }
    cat("\nCoefficients\n ~",
        if (is.na(p$dist)) {
          "flat"
        } else if (p$dist %in% c("normal", "student_t", "cauchy")) {
          if (is.null(p$df)) {
            paste0(p$dist, "(location = ", .fr2(p$location), 
                   ", scale = ", .fr2(p$scale), ")")
          } else {
            paste0(p$dist, "(df = ", .fr2(p$df), 
                  ", location = ", .fr2(p$location), 
                   ", scale = ", .fr2(p$scale),")")
          }
        } else if (p$dist %in% c("hs_plus")) {
          paste0("hs_plus(df1 = ", .fr2(p$df1), ", df2 = ", .fr2(p$df2), ")")
        } else if (p$dist %in% c("hs")) {
          paste0("hs(df = ", .fr2(p$df), ")")
        } else if (p$dist %in% c("R2")) {
          paste0("R2(location = ", .fr2(p$location), ", what = '", p$what, "')")
    })
    
    if (!is.null(p$adjusted_scale))
      cat("\n     **adjusted scale =", .fr3(p$adjusted_scale))
  }
  
  if (!is.null(prior_covariance)) {
    p <- prior_covariance
    p$regularization <- .format_pars(p$regularization, .fr2)
    p$concentration <- .format_pars(p$concentration, .fr2)
    p$shape <- .format_pars(prior_covariance$shape, .fr2)
    p$scale <- .format_pars(p$scale, .fr2)
    cat("\nCovariance\n ~",
        paste0(p$dist, "(",  "reg = ", .fr2(p$regularization), 
               ", conc = ", .fr2(p$concentration), ", shape = ", .fr2(p$shape), 
               ", scale = ", .fr2(p$scale), ")")
      )
  }
  
  if (!is.null(prior_counts)) {# stan_polr
    p <- prior_counts
    p$concentration <- .format_pars(p$concentration, .fr2)
    cat("\nCounts\n ~",
        paste0(p$dist, "(", "concentration = ", .fr2(p$concentration), ")"))
  }
  
  if (!is.null(x$scobit_exponent)) {
    p <- x$scobit_exponent
    cat("\nScobit Exponent\n ~",
        paste0(p$dist, "(shape = ", .fr2(p$shape), 
               ", rate = ", .fr2(p$rate), ")"))
  }
  
  invisible(x)
}


# internal ----------------------------------------------------------------
#
# @param x numeric vector
# @param formatter a formatting function to apply (see .fr2, .fr3 above)
# @param N the maximum number of values to include before replacing the rest
#   with '...'
.format_pars <- function(x, formatter, N = 3) {
  K <- length(x)
  if (K < 2)
    return(x)
  paste0(
    "[", 
    paste(c(formatter(x[1:min(N, K)]), if (N < K) "..."), 
          collapse = ","), 
    "]"
  )
}
