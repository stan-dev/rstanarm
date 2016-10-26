#' Extract and/or print a summary of the priors used for an rstanarm model
#'
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' 
#' @details 
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
prior_summary <- function(object, ...) {
  UseMethod("prior_summary")
}

#' @rdname prior_summary
#' @export
prior_summary.stanreg <- function(object, ...) {
  x <- object[["prior.info"]]
  if (is.null(x)) {
    message("Priors not found in stanreg object.")
    return(NULL)
  }
  structure(x, class = "prior_summary.stanreg", 
            model_name = deparse(substitute(object)))
}

#' @export
print.prior_summary.stanreg <- function(x, digits = 2, ...) {
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
    int_dist <- prior_intercept$dist
    cat("\nIntercept:\n ",
        if (is.na(int_dist)) {
          "flat"
        } else if (is.null(prior_intercept$df)) {
          with(prior_intercept, 
               paste0(dist,"(location = ", .fr2(location), 
                      ", scale = ", .fr2(scale),")"))
        } else {
          with(prior_intercept,
            paste0(dist, "(df = ", .fr2(df), ", location = ", .fr2(location), 
                   ", scale = ", .fr2(scale), ")"))
        }
      )
    if (!is.null(prior_intercept$adjusted_scale))
      cat("\n   adjusted scale =", .fr3(prior_intercept$adjusted_scale))
  }
  
  if (!is.null(prior_coef)) {
    coef_dist <- prior_coef$dist
    k <- length(prior_coef$location)
    if (!(coef_dist %in% c("R2", NA)) && k >= 2) {
      if (coef_dist %in% c("normal", "student_t", "cauchy")) {
        prior_coef$location <- .format_pars(prior_coef$location, .fr2)
        prior_coef$scale <- .format_pars(prior_coef$scale, .fr2)
        if (!is.null(prior_coef$df))
          prior_coef$df <- .format_pars(prior_coef$df, .fr2)
        if (!is.null(prior_coef$adjusted_scale))
          prior_coef$adjusted_scale <- .format_pars(prior_coef$adjusted_scale, .fr2)
      } else if (coef_dist %in% c("hs_plus")) {
        prior_coef$df1 <- .format_pars(prior_coef$df, .fr2)
        prior_coef$df2 <- .format_pars(prior_coef$scale, .fr2)
      } else if (coef_dist %in% c("hs")) {
        prior_coef$df <- .format_pars(prior_coef$df, .fr2)
      }
    }
    cat("\nCoefficients:\n ",
        if (is.na(coef_dist)) {
          "flat"
        } else if (coef_dist %in% c("normal", "student_t", "cauchy")) {
          if (is.null(prior_coef$df)) {
            with(prior_coef, 
                 paste0(dist, "(location = ", .fr2(location), 
                        ", scale = ", .fr2(scale), ")"))
          } else {
            with(prior_coef, 
                 paste0(dist, "(df = ", .fr2(df), ", location = ", .fr2(location), 
                        ", scale = ", .fr2(scale),")"))
          }
        } else if (coef_dist %in% c("hs_plus")) {
          with(prior_coef, 
               paste0("hs_plus(df1 = ", .fr2(df1), ", df2 = ", .fr2(df2), ")"))
        } else if (coef_dist %in% c("hs")) {
          with(prior_coef, 
               paste0("hs(df = ", .fr2(df), ")"))
        } else if (coef_dist %in% c("R2")) {
          with(prior_coef, 
               paste0("R2(location = ", .fr2(location), ", what = '", what, "')"))
    })
    
    if (!is.null(prior_coef$adjusted_scale))
      cat("\n   adjusted scale =", .fr3(prior_coef$adjusted_scale))
  }
  
  
  if (!is.null(prior_covariance)) {
    j <- length(prior_covariance$regularization)
    if (j >= 2) {
      prior_covariance$regularization <- .format_pars(prior_covariance$regularization, .fr2)
      prior_covariance$concentration <- .format_pars(prior_covariance$concentration, .fr2)
      prior_covariance$shape <- .format_pars(prior_covariance$shape, .fr2)
      prior_covariance$scale <- .format_pars(prior_covariance$scale, .fr2)
    }
    cat("\n Covariance:\n ",
        with(prior_covariance,
          paste0(dist, "(",  "reg = ", .fr2(regularization), 
                 ", conc = ", .fr2(concentration), ", shape = ", .fr2(shape), 
                 ", scale = ", .fr2(scale), ")")
        )
      )
  }
  
  if (!is.null(prior_counts)) {
    # stan_polr
    j <- length(prior_counts$concentration)
    if (j >= 2) {
      prior_counts$concentration <- .format_pars(prior_counts$concentration, .fr2)
    }
    cat("\n Counts:\n ",
        with(prior_counts, paste0(dist, "(", "conc = ", .fr2(concentration), ")"))
      )
  }
  
  invisible(x)
}


.format_pars <- function(p, formatter) {
  paste0(
    "[", 
    paste(c(formatter(p[1:2]), if (length(p) > 2) "..."), 
          collapse = ","), 
    "]"
  )
}
