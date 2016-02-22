#' @param prior The prior distribution for the regression coefficients. 
#'   \code{prior} can be a call to \code{normal}, \code{student_t},
#'   \code{cauchy}, \code{hs} or \code{hs_plus}. See \code{\link{priors}} for
#'   details. To to omit a prior ---i.e., to use a flat (improper) uniform
#'   prior--- set \code{prior} to \code{NULL}.
#' @param prior_intercept The prior distribution for the intercept.
#'   \code{prior_intercept} can be a call to \code{normal}, \code{student_t} or
#'   \code{cauchy}. See \code{\link{priors}} for details. To to omit a prior
#'   ---i.e., to use a flat (improper) uniform prior--- set
#'   \code{prior_intercept} to \code{NULL}. (\strong{Note:} the prior
#'   distribution for the intercept is set so it applies to the value when all
#'   predictors are centered.)
#' @param prior_ops Additional options related to prior distributions. Can 
#'   be \code{NULL} to omit a prior on the dispersion and see 
#'   \code{\link{prior_options}} otherwise.
