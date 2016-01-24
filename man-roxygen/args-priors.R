#' @param prior Prior for coefficients. Can be \code{NULL} to omit a prior and 
#'   see \code{\link{priors}} otherwise.
#' @param prior_intercept Prior for intercept. Can be \code{NULL} to omit a 
#'   prior and see \code{\link{priors}} otherwise. (\strong{Note:} the prior
#'   distribution for the intercept is set so it applies to the value when all
#'   predictors are centered.)
#' @param prior_ops Additional options related to prior distributions, specified
#'   via a call to \code{\link{prior_options}}. Among other things, 
#'   \code{prior_ops} can be used to set the prior scale on the dispersion in
#'   \code{\link{gaussian}} models and the shape and rate parameters of a gamma
#'   prior on the degrees of freedom when \code{family} is 
#'   \code{\link{t_family}}.
#'  
