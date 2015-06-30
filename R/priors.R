#' Prior distributions
#' 
#' These functions are used to specify the \code{prior} and 
#' \code{prior.for.intercept} arguments of the \code{stan_regression}, 
#' \code{stan_lm}, and \code{stan_glm} functions. 
#' 
#' @export 
#' @name priors
#' @param location prior location. Defaults to 0. For \code{normal} and 
#' \code{student_t} (provided that \code{df > 1}) this is the prior mean. For 
#' \code{cauchy} the mean is undefined and \code{location} is the prior median.
#' @param scale prior scale. Default depends (see Details). 
#' @param df prior degrees of freedom. Defaults to 1, in which case
#' \code{student_t} is equivalent to \code{cauchy}. 
#' @details If \code{scale} is not specified it will default to 10 for the 
#' intercept and 2.5 for the other coefficients, unless the probit link function 
#' is used, in which case these defaults are scaled by a factor of 
#' \code{dnorm(0)/dlogis(0)} (roughly 1.6).
#'
#' 
#' @examples
#' \dontrun{
#' stan_lm(y ~ x1 + x2, prior = student_t(4, 0, 2.5), prior.for.intercept = cauchy(0,10))
#' }
#' 
normal <- function(location = 0, scale = NULL) {
  validate_loc_scale_df(location, scale)
  nlist(dist = "normal", df = NA, location, scale)
}

#' @rdname priors
#' @export
student_t <- function(df = 1, location = 0, scale = NULL) {
  validate_loc_scale_df(location, scale, df)
  nlist(dist = "t", df, location, scale)
}

#' @rdname priors
#' @export
cauchy <- function(location = 0, scale = NULL) {
  student_t(df = 1, location = location, scale = scale)
}


