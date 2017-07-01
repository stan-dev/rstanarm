#' @param prior_intercept The prior distribution for the intercept. 
#'   \code{prior_intercept} can be a call to \code{normal}, \code{student_t} or 
#'   \code{cauchy}. See the \link[=priors]{priors help page} for details on 
#'   these functions. To omit a prior on the intercept ---i.e., to use a flat
#'   (improper) uniform prior--- \code{prior_intercept} can be set to
#'   \code{NULL}.
#'   
#'   \strong{Note:} If using a dense representation of the design matrix 
#'   ---i.e., if the \code{sparse} argument is left at its default value of
#'   \code{FALSE}--- then the prior distribution for the intercept is set so it
#'   applies to the value when all predictors are centered.
