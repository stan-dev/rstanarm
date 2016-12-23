#' @param prior The prior distribution for the regression coefficients. 
#'   \code{prior} can be a call to \code{normal}, \code{student_t}, 
#'   \code{cauchy}, \code{hs} or \code{hs_plus}. See \code{\link{priors}} for 
#'   details. To omit a prior ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior} can be set to \code{NULL}, although this is rarely a good 
#'   idea. (\strong{Note:} Unless \code{QR=TRUE}, if \code{prior} is specified 
#'   as \code{normal}, \code{student_t}, or \code{cauchy} with the 
#'   \code{autoscale} argument left at its default and recommended value of 
#'   \code{TRUE}, then the scale(s) of \code{prior} may be tuned internally 
#'   based on the scales of the predictors. See \code{\link{priors}} for details
#'   on the rescaling and \code{\link{prior_summary}} for a summary of the 
#'   priors used for a particular model.)
#' @param prior_intercept The prior distribution for the intercept. 
#'   \code{prior_intercept} can be a call to \code{normal}, \code{student_t} or 
#'   \code{cauchy}. See \code{\link{priors}} for details. To to omit a prior 
#'   ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_intercept} to \code{NULL}. (\strong{Note:} if a dense 
#'   representation of the design matrix is utilized ---i.e., if the 
#'   \code{sparse} argument is left at its default value of \code{FALSE}--- then
#'   the prior distribution for the intercept is set so it applies to the value 
#'   when all predictors are centered.)
#' @param prior_nuisance The prior distribution for the "nuisance" parameter (if
#'   applicable). The "nuisance" parameter refers to a different parameter 
#'   depending on the \code{family}. For Gaussian models it refers to the error 
#'   standard deviation \code{"sigma"}. For negative binomial models it is the 
#'   so-called \code{"size"} parameter (see e.g., \code{\link[stats]{rnbinom}}),
#'   and smaller values of \code{"size"} correspond to greater (over)dispersion.
#'   For gamma models it is the \code{"shape"} parameter (see e.g., 
#'   \code{\link[stats]{rgamma}}), and for inverse-Gaussian models it is the 
#'   so-called \code{"lambda"} parameter (which is essentially the reciprocal of
#'   a scale parameter). Binomial and Poisson models do not have nuisance 
#'   parameters. \code{prior_nuisance} can be a call to \code{exponential} to 
#'   use an exponential distribution, or \code{normal}, \code{student_t} or 
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy 
#'   prior. See \code{\link{priors}} for details on these functions. To omit a 
#'   prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_nuisance} to \code{NULL}.
