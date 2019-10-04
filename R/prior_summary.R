#' Summarize the priors used for an rstanarm model
#' 
#' The \code{prior_summary} method provides a summary of the prior distributions
#' used for the parameters in a given model. In some cases the user-specified
#' prior does not correspond exactly to the prior used internally by
#' \pkg{rstanarm} (see the sections below). Especially in these cases, but also
#' in general, it can be much more useful to visualize the priors. Visualizing
#' the priors can be done using the \code{\link{posterior_vs_prior}} function,
#' or alternatively by fitting the model with the \code{prior_PD} argument set
#' to \code{TRUE} (to draw from the prior predictive distribution instead of
#' conditioning on the outcome) and then plotting the parameters.
#' 
#' @aliases prior_summary
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param digits Number of digits to use for rounding.
#' @param ... Currently ignored by the method for stanreg objects.
#' 
#' @section Intercept (after predictors centered): 
#'   For \pkg{rstanarm} modeling functions that accept a \code{prior_intercept} 
#'   argument, the specified prior for the intercept term applies to the 
#'   intercept after \pkg{rstanarm} internally centers the predictors so they 
#'   each have mean zero. The estimate of the intercept returned to the user 
#'   correspond to the intercept with the predictors as specified by the user 
#'   (unmodified by \pkg{rstanarm}), but when \emph{specifying} the prior the 
#'   intercept can be thought of as the expected outcome when the predictors are
#'   set to their means. The only exception to this is for models fit with the 
#'   \code{sparse} argument set to \code{TRUE} (which is only possible with a
#'   subset of the modeling functions and never the default).
#'   
#' @section Adjusted scales: For some models you may see "\code{adjusted scale}"
#'   in the printed output and adjusted scales included in the object returned 
#'   by \code{prior_summary}. These adjusted scale values are the prior scales 
#'   actually used by \pkg{rstanarm} and are computed by adjusting the prior 
#'   scales specified by the user to account for the scales of the predictors 
#'   (as described in the documentation for the \code{\link[=priors]{autoscale}}
#'   argument). To disable internal prior scale adjustments set the 
#'   \code{autoscale} argument to \code{FALSE} when setting a prior using one of
#'   the distributions that accepts an \code{autoscale} argument. For example,
#'   \code{normal(0, 5, autoscale=FALSE)} instead of just \code{normal(0, 5)}.
#' 
#' @section Coefficients in Q-space:
#'   For the models fit with an \pkg{rstanarm} modeling function that supports 
#'   the \code{QR} argument (see e.g, \code{\link{stan_glm}}), if \code{QR} is 
#'   set to \code{TRUE} then the prior distributions for the regression
#'   coefficients specified using the \code{prior} argument are not relative to
#'   the original predictor variables \eqn{X} but rather to the variables in the
#'   matrix \eqn{Q} obtained from the \eqn{QR} decomposition of \eqn{X}. 
#'   
#'   In particular, if \code{prior = normal(location,scale)}, then this prior on
#'   the coefficients in \eqn{Q}-space can be easily translated into a joint 
#'   multivariate normal (MVN) prior on the coefficients on the original 
#'   predictors in \eqn{X}. Letting \eqn{\theta} denote the coefficients on
#'   \eqn{Q} and \eqn{\beta} the coefficients on \eqn{X} then if \eqn{\theta
#'   \sim N(\mu, \sigma)}{\theta ~ N(\mu, \sigma)} the corresponding prior on
#'   \eqn{\beta} is \eqn{\beta \sim MVN(R\mu, R'R\sigma^2)}{\beta ~ MVN(R\mu,
#'   R'R\sigma)}, where \eqn{\mu} and \eqn{\sigma} are vectors of the
#'   appropriate length. Technically, \pkg{rstanarm} uses a scaled \eqn{QR}
#'   decomposition to ensure that the columns of the predictor matrix used to
#'   fit the model all have unit scale, when the \code{autoscale} argument
#'   to the function passed to the \code{prior} argument is \code{TRUE} (the
#'   default), in which case the matrices actually used are
#'   \eqn{Q^\ast = Q \sqrt{n-1}}{Q* = Q (n-1)^0.5} and \eqn{R^\ast =
#'   \frac{1}{\sqrt{n-1}} R}{R* = (n-1)^(-0.5) R}. If \code{autoscale = FALSE}
#'   we instead scale such that the lower-right element of \eqn{R^\ast}{R*} is 
#'   \eqn{1}, which is useful if you want to specify a prior on the coefficient 
#'   of the last predictor in its original units (see the documentation for the 
#'   \code{\link[=stan_glm]{QR}} argument).
#'   
#'   If you are interested in the prior on \eqn{\beta} implied by the prior on
#'   \eqn{\theta}, we strongly recommend visualizing it as described above in
#'   the \strong{Description} section, which is simpler than working it out
#'   analytically.
#'   
#' @return A list of class "prior_summary.stanreg", which has its own print
#'   method.
#'   
#' @seealso The \link[=priors]{priors help page} and the \emph{Prior
#'   Distributions} vignette.
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
#' # autoscale=FALSE for prior on coefficients
#' fit <- stan_glm(mpg ~ wt + am, data = mtcars, 
#'                 prior = normal(0, c(2.5, 4)), 
#'                 prior_intercept = normal(0, 5), 
#'                 iter = 10, chains = 1) # only for demonstration 
#' prior_summary(fit)
#' 
#' fit2 <- update(fit, prior = normal(0, c(2.5, 4), autoscale=FALSE), 
#'                prior_intercept = normal(0, 5, autoscale=FALSE))
#' prior_summary(fit2)
#' 
prior_summary.stanreg <- function(object, digits = 2,...) {
  x <- object[["prior.info"]]
  if (is.null(x)) {
    message("Priors not found in stanreg object.")
    return(invisible(NULL))
  }  
  if (is.stanmvreg(object)) {
    M <- get_M(object)
    x <- structure(x, M = M) 
  }
  structure(x, class = "prior_summary.stanreg",
            QR = used.QR(object),
            sparse = used.sparse(object),
            model_name = deparse(substitute(object)), 
            stan_function = object$stan_function,
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
  QR <- attr(x, "QR")
  sparse <- attr(x, "sparse")
  model_name <- attr(x, "model_name")
  stan_function <- attr(x, "stan_function")
  
  msg <- paste0("Priors for model '", model_name, "'")
  cat(msg, "\n------")
  
  if (!stan_function == "stan_mvmer") {
    if (!is.null(x[["prior_intercept"]]))
      .print_scalar_prior(
        x[["prior_intercept"]], 
        txt = paste0("Intercept", if (!sparse) " (after predictors centered)"), 
        formatters
      )
    if (!is.null(x[["prior"]]))
      .print_vector_prior(
        x[["prior"]], 
        txt = paste0("\nCoefficients", if (QR) " (in Q-space)"), 
        formatters = formatters
      )
    if (!is.null(x[["prior_aux"]])) {
      aux_name <- x[["prior_aux"]][["aux_name"]]
      aux_dist <- x[["prior_aux"]][["dist"]]
      if (aux_dist %in% c("normal", "student_t", "cauchy"))
        x[["prior_aux"]][["dist"]] <- paste0("half-", aux_dist)
      .print_scalar_prior(
        x[["prior_aux"]], 
        txt = paste0("\nAuxiliary (", aux_name, ")"), 
        formatters
      )
    }    
  } else { # unique to stan_mvmer
    M <- attr(x, "M")
    for (m in 1:M) {
      if (!is.null(x[["prior_intercept"]][[m]]))
        .print_scalar_prior(
          x[["prior_intercept"]][[m]], 
          txt = paste0(if (m > 1) "\n", "y", m, "|Intercept", if (!sparse) 
            " (after predictors centered)"), 
          formatters
        )
      if (!is.null(x[["prior"]][[m]]))
        .print_vector_prior(
          x[["prior"]][[m]], 
          txt = paste0("\ny", m, "|Coefficients", if (QR) " (in Q-space)"), 
          formatters = formatters
        )
      if (!is.null(x[["prior_aux"]][[m]])) {
        aux_name <- x[["prior_aux"]][[m]][["aux_name"]]
        aux_dist <- x[["prior_aux"]][[m]][["dist"]]
        if (aux_dist %in% c("normal", "student_t", "cauchy"))
          x[["prior_aux"]][[m]][["dist"]] <- paste0("half-", aux_dist)
        .print_scalar_prior(
          x[["prior_aux"]][[m]], 
          txt = paste0("\ny", m, "|Auxiliary (", aux_name, ")"), 
          formatters
        )
      }
    }    
  }
  
  # unique to stan_betareg
  if (!is.null(x[["prior_intercept_z"]]))
    .print_scalar_prior(
      x[["prior_intercept_z"]], 
      txt = paste0("\nIntercept_z", if (!sparse) " (after predictors centered)"), 
      formatters
    )
  if (!is.null(x[["prior_z"]]))
    .print_vector_prior(x[["prior_z"]], txt = "\nCoefficients_z", formatters)
  
  # unique to stan_jm
  if (stan_function == "stan_jm") {
    M <- attr(x, "M")
    for (m in 1:M) {
      if (!is.null(x[["priorLong_intercept"]][[m]]))
        .print_scalar_prior(
          x[["priorLong_intercept"]][[m]], 
          txt = paste0(if (m > 1) "\n", "Long", m, "|Intercept", if (!sparse) 
            " (after predictors centered)"), 
          formatters
        )
      if (!is.null(x[["priorLong"]][[m]]))
        .print_vector_prior(
          x[["priorLong"]][[m]], 
          txt = paste0("\nLong", m, "|Coefficients", if (QR) " (in Q-space)"), 
          formatters = formatters
        )
      if (!is.null(x[["priorLong_aux"]][[m]])) {
        aux_name <- x[["priorLong_aux"]][[m]][["aux_name"]]
        aux_dist <- x[["priorLong_aux"]][[m]][["dist"]]
        if (aux_dist %in% c("normal", "student_t", "cauchy"))
          x[["priorLong_aux"]][[m]][["dist"]] <- paste0("half-", aux_dist)
        .print_scalar_prior(
          x[["priorLong_aux"]][[m]], 
          txt = paste0("\nLong", m, "|Auxiliary (", aux_name, ")"), 
          formatters
        )
      }
    }
    if (!is.null(x[["priorEvent_intercept"]]))
      .print_scalar_prior(
        x[["priorEvent_intercept"]], 
        txt = paste0("\nEvent|Intercept", if (!sparse) " (after predictors centered)"), 
        formatters
      )
    if (!is.null(x[["priorEvent"]]))
      .print_vector_prior(
        x[["priorEvent"]], 
        txt = "\nEvent|Coefficients", 
        formatters = formatters
      )
    if (!is.null(x[["priorEvent_aux"]])) {
      aux_name <- x[["priorEvent_aux"]][["aux_name"]]
      aux_dist <- x[["priorEvent_aux"]][["dist"]]
      if ((aux_name == "weibull-shape") &&
          (aux_dist %in% c("normal", "student_t", "cauchy"))) { # weibull
        x[["priorEvent_aux"]][["dist"]] <- paste0("half-", aux_dist)
        .print_scalar_prior(
          x[["priorEvent_aux"]], 
          txt = paste0("\nEvent|Auxiliary (", aux_name, ")"), 
          formatters
        )        
      } else { # bs or piecewise
        .print_vector_prior(
          x[["priorEvent_aux"]], 
          txt = paste0("\nEvent|Auxiliary (", aux_name, ")"), 
          formatters
        )
      }
    }
    if (!is.null(x[["priorEvent_assoc"]]))
      .print_vector_prior(
        x[["priorEvent_assoc"]], 
        txt = "\nAssociation parameters", 
        formatters = formatters
      )
  }  
  
  # unique to stan_(g)lmer, stan_gamm4, stan_mvmer, or stan_jm
  if (!is.null(x[["prior_covariance"]]))
    .print_covariance_prior(x[["prior_covariance"]], txt = "\nCovariance", formatters)
  
  # unique to stan_polr
  if (!is.null(x[["prior_counts"]])) {
    p <- x[["prior_counts"]]
    p$concentration <- .format_pars(p$concentration, .fr2)
    cat("\n\nCounts\n ~",
        paste0(p$dist, "(", "concentration = ", .fr2(p$concentration), ")"))
  }
  if (!is.null(x[["scobit_exponent"]])) {
    p <- x[["scobit_exponent"]]
    cat("\n\nScobit Exponent\n ~",
        paste0(p$dist, "(shape = ", .fr2(p$shape), 
               ", rate = ", .fr2(p$rate), ")"))
  }

  cat("\n------\n")
  cat("See help('prior_summary.stanreg') for more details\n")
  invisible(x)
}


# internal ----------------------------------------------------------------

# check if model was fit using QR=TRUE
used.QR <- function(x) {
  isTRUE(getCall(x)[["QR"]])
}

# check if model was fit using sparse=TRUE
used.sparse <- function(x) {
  isTRUE(getCall(x)[["sparse"]])
}

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
  
  .cat_scalar_prior <- function(p, adjusted = FALSE, prepend_chars = "\n ~") {
    if (adjusted) {
      p$scale <- p$adjusted_scale
      p$rate <- 1/p$adjusted_scale
    }
    cat(prepend_chars,
        if (is.na(p$dist)) {
          "flat"
        } else if (p$dist == "exponential") {
          paste0(p$dist,"(rate = ", .f1(p$rate), ")")
        } else { # normal, student_t, cauchy
          if (is.null(p$df)) {
            paste0(p$dist,"(location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale),")")
          } else {
            paste0(p$dist, "(df = ", .f1(p$df), 
                   ", location = ", .f1(p$location), 
                   ", scale = ", .f1(p$scale), ")")
          }
        }
    )
  }
  cat(paste0("\n", txt))
  if (is.null(p$adjusted_scale)) {
    .cat_scalar_prior(p, adjusted = FALSE)
  } else {
    cat("\n  Specified prior:")
    .cat_scalar_prior(p, adjusted = FALSE, prepend_chars = "\n    ~")
    cat("\n  Adjusted prior:")
    .cat_scalar_prior(p, adjusted = TRUE, prepend_chars =  "\n    ~")
  }
  
}

.print_covariance_prior <- function(p, txt = "Covariance", formatters = list()) {
  if (p$dist == "decov") {
    .f1 <- formatters[[1]]
    p$regularization <- .format_pars(p$regularization, .f1)
    p$concentration <- .format_pars(p$concentration, .f1)
    p$shape <- .format_pars(p$shape, .f1)
    p$scale <- .format_pars(p$scale, .f1)
    cat(paste0("\n", txt, "\n ~"),
        paste0(p$dist, "(",  
               "reg. = ",    .f1(p$regularization),
               ", conc. = ", .f1(p$concentration), 
               ", shape = ", .f1(p$shape),
               ", scale = ", .f1(p$scale), ")")
    )    
  } else if (p$dist == "lkj") {
    .f1 <- formatters[[1]]
    .f2 <- formatters[[2]]
    p$regularization <- .format_pars(p$regularization, .f1)
    p$df <- .format_pars(p$df, .f1)
    p$scale <- .format_pars(p$scale, .f1)
    if (!is.null(p$adjusted_scale))
      p$adjusted_scale <- .format_pars(p$adjusted_scale, .f2)
    cat(paste0("\n", txt, "\n ~"),
        paste0(p$dist, "(",  
               "reg. = ",    .f1(p$regularization),
               ", df = ",    .f1(p$df), 
               ", scale = ", .f1(p$scale), ")")
    )    
    if (!is.null(p$adjusted_scale))
      cat("\n     **adjusted scale =", .f2(p$adjusted_scale))
  }
}



.print_vector_prior <- function(p, txt = "Coefficients", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  
  if (!(p$dist %in% c("R2", NA))) {
    if (p$dist %in% c("normal", "student_t", "cauchy", "laplace", "lasso", "product_normal")) {
      p$location <- .format_pars(p$location, .f1)
      p$scale <- .format_pars(p$scale, .f1)
      if (!is.null(p$df))
        p$df <- .format_pars(p$df, .f1)
      if (!is.null(p$adjusted_scale))
        p$adjusted_scale <- .format_pars(p$adjusted_scale, .f2)
    } else if (p$dist %in% c("hs_plus")) {
      p$df1 <- .format_pars(p$df, .f1)
      p$df2 <- .format_pars(p$scale, .f1)
    } else if (p$dist %in% c("hs")) {
      p$df <- .format_pars(p$df, .f1)
    } else if (p$dist %in% c("product_normal"))
      p$df <- .format_pars(p$df, .f1)
  }
  
  .cat_vector_prior <- function(p, adjusted = FALSE, prepend_chars = "\n ~") {
    if (adjusted) {
      p$scale <- p$adjusted_scale
    }
    cat(prepend_chars, 
        if (is.na(p$dist)) {
          "flat"
        } else if (p$dist %in% c("normal", "student_t", "cauchy", 
                                 "laplace", "lasso", "product_normal")) {
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
  }
  
  cat(paste0("\n", txt))
  if (is.null(p$adjusted_scale)) {
    .cat_vector_prior(p, adjusted = FALSE)
  } else {
    cat("\n  Specified prior:")
    .cat_vector_prior(p, adjusted = FALSE, prepend_chars = "\n    ~")
    cat("\n  Adjusted prior:")
    .cat_vector_prior(p, adjusted = TRUE, prepend_chars =  "\n    ~")
  }
}
