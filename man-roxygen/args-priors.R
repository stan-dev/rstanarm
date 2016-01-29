#' @param prior Prior for coefficients. See \code{\link{priors}} for details. 
#'   Set \code{prior} to \code{NULL} to omit a prior, i.e., use an (improper)
#'   uniform prior.
#' @param prior_intercept Prior for intercept. See \code{\link{priors}} for 
#'   details. Set \code{prior_intercept} to \code{NULL} to omit a prior, i.e.,
#'   use an (improper) uniform prior. (\strong{Note:} the prior distribution for
#'   the intercept is set so it applies to the value when all predictors are 
#'   centered.)
