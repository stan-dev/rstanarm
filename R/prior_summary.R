#' Extract and/or print a summary of the priors used for an rstanarm model
#'
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' 
#' @return A list of class "prior_summary.stanreg", which has its own print
#'   method.
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
  structure(x, class = "prior_summary.stanreg")
}

#' @export
print.prior_summary.stanreg <- function(x, digits = 2, ...) {
  priors <- x
  .dig <- digits
  .fr2 <- function(x, .digits=.dig, ...) format(x, digits = .digits, ...)
  .fr3 <- function(x, .nsmall=.dig) .fr2(x, nsmall = .nsmall)
  
  prior_intercept <- priors[["prior_intercept"]]
  prior_coef <- priors[["prior"]]
  prior_covariance <- priors[["prior_covariance"]]
  prior_counts <- priors[["prior_counts"]] # for stan_polr
  
  cat("Priors:")
  if (!is.null(prior_intercept)) {
    int_dist <- prior_intercept$dist
    cat("\n Intercept:",
        if (is.na(int_dist)) {
          "flat"
        } else if (is.null(prior_intercept$df)) {
          with(prior_intercept,
               paste0(dist, "(loc = ", .fr2(location), ", scale = ", .fr2(scale), ")"))
        } else {
          with(
            prior_intercept,
            paste0(dist, "(df = ", .fr2(df), ", loc = ", .fr2(location), ", scale = ", .fr2(scale), ")")
          )
        })
    if (!is.null(prior_intercept$adjusted_scale))
      cat("\n  adjusted scale =", .fr3(prior_intercept$adjusted_scale))
  }
  
  if (!is.null(prior_coef)) {
    coef_dist <- prior_coef$dist 
    k <- length(prior_coef$location)
    if (!(coef_dist %in% c("R2", NA)) && k >= 2) {
      if (coef_dist %in% c("normal", "student_t", "cauchy")) {
        if (!is.null(prior_coef$df))
          prior_coef$df <-
            paste0("[", paste(c(.fr2(prior_coef$df[1:2]), if (k > 2)
              "..."), collapse = ","), "]")
        prior_coef$location <-
          paste0("[", paste(c(.fr2(prior_coef$location[1:2]), if (k > 2)
            "..."), collapse = ","), "]")
        prior_coef$scale <-
          paste0("[", paste(c((prior_coef$scale[1:2]), if (k > 2)
            "..."), collapse = ","), "]")
        if (!is.null(prior_coef$adjusted_scale))
          prior_coef$adjusted_scale <-
          paste0("[", paste(c(.fr2(prior_coef$adjusted_scale[1:2]), if (k > 2)
            "..."), collapse = ","), "]")
      } else if (coef_dist %in% c("hs_plus")) {
        prior_coef$df1 <-
          paste0("[", paste(c(.fr2(prior_coef$df[1:2]), if (k > 2)
            "..."), collapse = ","), "]")
        prior_coef$df2 <-
          paste0("[", paste(c(.fr2(prior_coef$scale[1:2]), if (k > 2)
            "..."), collapse = ","), "]")
      } else if (coef_dist %in% c("hs")) {
        prior_coef$df <-
          paste0("[", paste(c(.fr2(prior_coef$df[1:2]), if (k > 2)
            "..."), collapse = ","), "]")
      }
    }
    
    cat("\n Coefficients:",
        if (is.na(coef_dist)) {
          "flat"
        } else if (coef_dist %in% c("normal", "student_t", "cauchy")) {
          if (is.null(prior_coef$df)) {
            with(prior_coef, paste0(dist, "(loc = ", .fr2(location), ", scale = ", .fr2(scale), ")"))
          } else {
            with(prior_coef, paste0(dist, "(df = ", .fr2(df), ", loc = ", .fr2(location), ", scale = ", .fr2(scale), ")"))
          }
        } else if (coef_dist %in% c("hs_plus")) {
          with(prior_coef, paste0("hs_plus(df1 = ", .fr2(df1), ", df2 = ", .fr2(df2), ")"))
        } else if (coef_dist %in% c("hs")) {
          with(prior_coef, paste0("hs(df = ", .fr2(df), ")"))
        } else if (coef_dist %in% c("R2")) {
          with(prior_coef, paste0("R2(loc = ", .fr2(location), ", what = '", what, "')"))
        }
    )
    
    if (!is.null(prior_coef$adjusted_scale))
      cat("\n  adjusted scale(s) =", .fr3(prior_coef$adjusted_scale))
  }
  
  
  if (!is.null(prior_covariance)) {
    j <- length(prior_covariance$regularization)
    if (j >= 2) {
      prior_covariance$regularization <-
        paste0("[", paste(
          c(.fr2(prior_covariance$regularization[1:2]), if (j > 2)
            "..."),
          collapse = ","
        ), "]")
      prior_covariance$concentration <-
        paste0("[", paste(
          c(.fr2(prior_covariance$concentration[1:2]), if (j > 2)
            "..."),
          collapse = ","
        ), "]")
      prior_covariance$shape <-
        paste0("[", paste(c(
          .fr2(prior_covariance$shape[1:2]), if (j > 2)
            "..."
        ), collapse = ","), "]")
      prior_covariance$scale <-
        paste0("[", paste(c(
          .fr2(prior_covariance$scale[1:2]), if (j > 2)
            "..."
        ), collapse = ","), "]")
    }
    cat("\n Covariance:",
        with(
          prior_covariance,
          paste0(
            dist,
            "(",
            "reg = ",
            .fr2(regularization),
            ", conc = ",
            .fr2(concentration),
            ", shape = ",
            .fr2(shape),
            ", scale = ",
            .fr2(scale),
            ")"
          )
        ))
  }
  
  if (!is.null(prior_counts)) { # stan_polr
    j <- length(prior_counts$concentration)
    if (j >= 2) {
      prior_counts$concentration <-
        paste0("[", paste(
          c(.fr2(prior_counts$concentration[1:2]), if (j > 2)
            "..."),
          collapse = ","
        ), "]")
    }
    cat("\n Counts:",
        with(
          prior_counts,
          paste0(
            dist,
            "(",
            "conc = ",
            .fr2(concentration),
            ")"
          )
        ))
  }
}
