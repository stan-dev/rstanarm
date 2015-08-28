#' @param prior Prior for coefficients. Can be \code{NULL} to omit a prior and 
#'   see \code{\link{priors}} otherwise.
#' @param prior.for.intercept Prior for intercept. Can be \code{NULL} to omit a 
#'   prior and see \code{\link{priors}} otherwise. (\strong{Note:} the prior
#'   distribution for the intercept is set so it applies to the value when all
#'   predictors are centered.)
#' @param prior.options Additional options related to prior distributions. Can 
#'   be \code{NULL} to omit a prior on the dispersion and see 
#'   \code{\link{priors}} otherwise.
