#' Prior distributions
#' 
#' These functions are used to specify the \code{prior} and 
#' \code{prior.for.intercept}, and \code{prior.options} arguments of the
#' \code{stan_lm}, and \code{stan_glm} functions.
#' 
#' @export 
#' @name priors
#' @param location Prior location. Defaults to 0. For \code{normal} and 
#'   \code{student_t} (provided that \code{df > 1}) this is the prior mean. For 
#'   \code{cauchy} (which is equivalent to \code{student_t} with \code{df=1}), 
#'   the mean is undefined and \code{location} is the prior median.
#' @param scale Prior scale. Default depends (see Details).
#' @param df Prior degrees of freedom. Defaults to 1, in which case 
#'   \code{student_t} is equivalent to \code{cauchy}.
#'   
#' @details For the prior distribution for the intercept, \code{location}, 
#'   \code{scale}, and \code{df} should be scalars. For the prior for the other
#'   coefficients they can either be vectors of length equal to the number of
#'   coefficients (not including the intercept), or they can be scalars, in 
#'   which case they will be replicated to the appropriate length.
#'   
#'   If \code{scale} is not specified it will default to 10 for the intercept
#'   and 2.5 for the other coefficients, unless the probit link function is
#'   used, in which case these defaults are scaled by a factor of 
#'   \code{dnorm(0)/dlogis(0)} (roughly 1.6).
#'
#' @return A named list.
#' @examples
#' \dontrun{
#' stan_lm(y ~ x1 + x2, prior = student_t(4, 0, 2.5), prior.for.intercept = cauchy(0,10))
#' }
#' 
normal <- function(location = 0, scale = NULL) {
  validate_parameter_value(scale)
  nlist(dist = "normal", df = NA, location, scale)
}

#' @rdname priors
#' @export
student_t <- function(df = 1, location = 0, scale = NULL) {
  validate_parameter_value(scale)
  validate_parameter_value(df)
  nlist(dist = "t", df, location, scale)
}

#' @rdname priors
#' @export
cauchy <- function(location = 0, scale = NULL) {
  student_t(df = 1, location = location, scale = scale)
}

#' @rdname priors
#' @export 
#' @param prior.scale.for.dispersion Prior scale for the standard error of the 
#'   regression in Gaussian models, which is given a half-Cauchy prior truncated
#'   at zero.
#' @param min.prior.scale Minimum prior scale for the intercept and 
#'   coefficients.
#' @param scaled Logical, defaulting to \code{TRUE}. If \code{TRUE} the 
#'   \code{prior.scale} is further scaled by the range of the predictor if the 
#'   predictor has exactly two unique values and scales prior.scale by twice the
#'   standard deviation of the predictor if it has more than two unique values.
#'
prior_options <- function(prior.scale.for.dispersion = 5, 
                          min.prior.scale = 1e-12, 
                          scaled = TRUE) {
  validate_parameter_value(prior.scale.for.dispersion)
  validate_parameter_value(min.prior.scale)
  nlist(scaled, min.prior.scale, prior.scale.for.dispersion)
}

