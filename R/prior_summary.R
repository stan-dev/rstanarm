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
    return(invisible(NULL))
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
  formatters <- list(.fr2, .fr3)
  msg <- paste0("Priors for model '", attr(x, "model_name"), "'")
  cat(msg, "\n------")
  
  if (!is.null(x[["prior_intercept"]]))
    .print_scalar_prior(x[["prior_intercept"]], txt = "Intercept", formatters)
  if (!is.null(x[["prior"]]))
    .print_vector_prior(x[["prior"]], txt = "Coefficients", formatters)
  if (!is.null(x[["prior_dispersion"]])) {
    dispersion_name <- x[["prior_dispersion"]][["dispersion_name"]]
    .print_scalar_prior(x[["prior_dispersion"]], txt = dispersion_name, formatters)
  }
  
  # unique to stan_betareg
  if (!is.null(x[["prior_intercept_z"]]))
    .print_scalar_prior(x[["prior_intercept_z"]], txt = "Intercept_z", formatters)
  if (!is.null(x[["prior_z"]]))
    .print_vector_prior(x[["prior_z"]], txt = "Coefficients_z", formatters)
  
  # unique to stan_(g)lmer or stan_gamm4
  if (!is.null(x[["prior_covariance"]]))
    .print_covariance_prior(x[["prior_covariance"]], txt = "Covariance", formatters)
  
  # unique to stan_polr
  if (!is.null(x[["prior_counts"]])) {
    p <- x[["prior_counts"]]
    p$concentration <- .format_pars(p$concentration, .fr2)
    cat("\nCounts\n ~",
        paste0(p$dist, "(", "concentration = ", .fr2(p$concentration), ")"))
  }
  if (!is.null(x[["scobit_exponent"]])) {
    p <- x[["scobit_exponent"]]
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

# Print priors for intercept/coefs (called internally by print.prior_summary.stanreg)
#
# @param p named list of prior stuff
# @param txt header to be printed
# @param formatters a list of two formatter functions like .fr2, .fr3 (defined
#   in prior_summary.stanreg). The first is used for format all numbers except
#   for adjusted scales, for which the second function is used. This is kind of
#   hacky and should be replaced at some point.
# 
.print_scalar_prior <- function(p, txt = "Intercept", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  cat(paste0("\n", txt, "\n ~"),
      if (is.na(p$dist)) {
        "flat"
      } else if (is.null(p$df)) {
        paste0(p$dist,"(location = ", .f1(p$location), 
               ", scale = ", .f1(p$scale),")")
      } else {
        paste0(p$dist, "(df = ", .f1(p$df), ", 
               location = ", .f1(p$location), 
               ", scale = ", .f1(p$scale), ")")
      }
  )
  if (!is.null(p$adjusted_scale))
    cat("\n     **adjusted scale =", .f2(p$adjusted_scale))
      }

.print_vector_prior <- function(p, txt = "Coefficients", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  
  if (!(p$dist %in% c("R2", NA))) {
    if (p$dist %in% c("normal", "student_t", "cauchy")) {
      p$location <- .format_pars(p$location, .f1)
      p$scale <- .format_pars(p$scale, .f1)
      if (!is.null(p$df))
        p$df <- .format_pars(p$df, .f1)
      if (!is.null(p$adjusted_scale))
        p$adjusted_scale <- .format_pars(p$adjusted_scale, .f1)
    } else if (p$dist %in% c("hs_plus")) {
      p$df1 <- .format_pars(p$df, .f1)
      p$df2 <- .format_pars(p$scale, .f1)
    } else if (p$dist %in% c("hs")) {
      p$df <- .format_pars(p$df, .f1)
    }
  }
  cat(paste0("\n", txt, "\n ~"),
      if (is.na(p$dist)) {
        "flat"
      } else if (p$dist %in% c("normal", "student_t", "cauchy")) {
        if (is.null(p$df)) {
          paste0(p$dist, "(location = ", .f1(p$location), 
                 ", scale = ", .f1(p$scale), ")")
        } else {
          paste0(p$dist, "(df = ", .f1(p$df), 
                 ", location = ", .f1(p$location), 
                 ", scale = ", .f1(p$scale),")")
        }
      } else if (p$dist %in% c("hs_plus")) {
        paste0("hs_plus(df1 = ", .f1(p$df1), ", df2 = ", .f1(p$df2), ")")
      } else if (p$dist %in% c("hs")) {
        paste0("hs(df = ", .f1(p$df), ")")
      } else if (p$dist %in% c("R2")) {
        paste0("R2(location = ", .f1(p$location), ", what = '", p$what, "')")
      })
  
  if (!is.null(p$adjusted_scale))
    cat("\n     **adjusted scale =", .f2(p$adjusted_scale))
}
.print_covariance_prior <- function(p, txt = "Covariance", formatters = list()) {
  .f1 <- formatters[[1]]
  p$regularization <- .format_pars(p$regularization, .f1)
  p$concentration <- .format_pars(p$concentration, .f1)
  p$shape <- .format_pars(p$shape, .f1)
  p$scale <- .format_pars(p$scale, .f1)
  cat(paste0("\n", txt, "\n ~"),
      paste0(p$dist, "(",  "reg = ", .f1(p$regularization),
             ", conc = ", .f1(p$concentration), ", shape = ", .f1(p$shape),
             ", scale = ", .f1(p$scale), ")")
  )
}
