#' Deprecated functions
#' 
#' These functions are deprecated and will be removed in a future release. The 
#' \strong{Arguments} section below provides details on how the functionality 
#' obtained via each of the arguments has been replaced.
#' 
#' @name rstanarm-deprecated
#' 
NULL

#' @rdname rstanarm-deprecated
#' @export 
#' @param prior_scale_for_dispersion,min_prior_scale,scaled Arguments to 
#'   deprecated \code{prior_options} function. The functionality provided 
#'   by the now deprecated \code{prior_options} function has been replaced 
#'   as follows: 
#'   \describe{
#'   \item{\code{prior_scale_for_dispersion}}{
#'    Instead of using the \code{prior_scale_for_dispersion} argument to 
#'    \code{prior_options}, priors for these parameters can now be 
#'    specified directly when calling \code{\link{stan_glm}} (or
#'    \code{\link{stan_glmer}}, etc.) using the new \code{prior_aux}
#'    argument.
#'   }
#'   \item{\code{scaled}}{
#'    Instead of setting \code{prior_options(scaled=FALSE)}, internal rescaling
#'    is now toggled using the new \code{autoscale} arguments to
#'    \code{\link{normal}}, \code{\link{student_t}}, and \code{\link{cauchy}} 
#'    (the other prior distributions do not support 'autoscale').
#'   }
#'   \item{\code{min_prior_scale}}{
#'    No replacement. \code{min_prior_scale} (the minimum possible scale
#'    parameter value that be used for priors) is now fixed to \code{1e-12}.
#'   }
#'   }
#' 
prior_options <- function(prior_scale_for_dispersion = 5, 
                          min_prior_scale = 1e-12, 
                          scaled = TRUE) {
  warning(
    "'prior_options' is deprecated and will be removed in a future release.",
    "\n* Priors for auxiliary parameters should now be set using",
    " the new 'prior_aux' argument when calling ",
    "'stan_glm', 'stan_glmer', etc.",
    "\n* Instead of setting 'prior_options(scaled=FALSE)',", 
    " internal rescaling is now toggled using the", 
    " new 'autoscale' argument to 'normal', 'student_t', or 'cauchy'", 
    " (the other prior distributions do not support 'autoscale').", 
    call. = FALSE
  )
  validate_parameter_value(prior_scale_for_dispersion)
  validate_parameter_value(min_prior_scale)
  out <- nlist(scaled, min_prior_scale, prior_scale_for_dispersion)
  structure(out, from_prior_options = TRUE)
}


# function used in stan_glm.fit to preserve backwards compatibility.
# should be removed when prior_options is officially removed
.support_deprecated_prior_options <- 
  function(prior, 
           prior_intercept, 
           prior_aux, 
           prior_ops) {
    if (!isTRUE(attr(prior_ops, "from_prior_options")))
      stop(
        "The 'prior_ops' argument must be a call to 'prior_options'. ",
        "But 'prior_options' is deprecated and will be removed in a future release. ", 
        "See help('rstanarm-deprecated') for details on the functionality ", 
        "that replaces 'prior_options'.", 
        call. = FALSE
      )
    
    po_disp_scale <- prior_ops[["prior_scale_for_dispersion"]]
    po_scaled <- prior_ops[["scaled"]]
    
    if (!is.null(prior_aux) && !is.null(po_disp_scale)) {
      if (po_disp_scale != prior_aux[["scale"]]) {
        warning(
          "Setting prior scale for aux to value specified in ",
          "'prior_options' rather than value specified in 'prior_aux'.",
          call. = FALSE
        )
        prior_aux[["scale"]] <- po_disp_scale
      }
    }
    if (!is.null(po_scaled) && identical(po_scaled, FALSE)) {
      if (isTRUE(prior$dist %in% c("normal", "t")))
        prior$autoscale <- FALSE
      if (!is.null(prior_intercept))
        prior_intercept$autoscale <- FALSE
    }
    
    nlist(prior, prior_intercept, prior_aux)
  }

