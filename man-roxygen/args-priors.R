#' @param prior The prior distribution for the regression coefficients. 
#'   \code{prior} can be a call to \code{normal}, \code{student_t}, 
#'   \code{cauchy}, \code{hs} or \code{hs_plus}. See \code{\link{priors}} for 
#'   details. To omit a prior ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior} can be set to \code{NULL}, although this is rarely a good 
#'   idea. (\strong{Note:} unless \code{QR=TRUE}, if the \code{scaled} argument
#'   to \code{\link{prior_options}} is left at its default and recommended value
#'   of \code{TRUE}, then the scale(s) of \code{prior} may be modified
#'   internally based on the scales of the predictors, as in the \pkg{arm}
#'   package. See \code{\link{priors}} for details on the rescaling and 
#'   \code{\link{prior_summary}} for a summary of the priors used for a
#'   particular model.)
#' @param prior_intercept The prior distribution for the intercept. 
#'   \code{prior_intercept} can be a call to \code{normal}, \code{student_t} or 
#'   \code{cauchy}. See \code{\link{priors}} for details. To to omit a prior 
#'   ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_intercept} to \code{NULL}. (\strong{Note:} if a dense 
#'   representation of the design matrix is utilized ---i.e., if the
#'   \code{sparse} argument is left at its default value of \code{FALSE}--- then
#'   the prior distribution for the intercept is set so it applies to the value
#'   when all predictors are centered.)
#' @param prior_dispersion The prior distribution for the "dispersion" parameter
#'   (if applicable). The "dispersion" parameter refers to a different parameter
#'   depending on the \code{family}. For Gaussian models it is the residual SD 
#'   sigma, for negative binomial models it is the overdispersion parameter, for
#'   gamma models it is the shape parameter, and for inverse-Gaussian models it 
#'   is the lambda parameter. Binomial and Poisson models do not have dispersion
#'   parameters. \code{prior_dispersion} can be a call to \code{normal},
#'   \code{student_t} or \code{cauchy}, which results in a half-normal, half-t,
#'   or half-Cauchy prior. See \code{\link{priors}} for details. To to omit a
#'   prior ---i.e., to use a flat (improper) uniform prior--- set
#'   \code{prior_dispersion} to \code{NULL}.
#' 
