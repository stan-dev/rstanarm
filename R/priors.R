#' Prior distributions
#' 
#' These functions are used to specify the \code{prior}, 
#' \code{prior.for.intercept}, and \code{prior.options} arguments of various
#' model-fitting functions in this package.
#' 
#' @export 
#' @name priors
#' @param location Prior location. For \code{normal} and 
#'   \code{student_t} (provided that \code{df > 1}) this is the prior mean. For 
#'   \code{cauchy} (which is equivalent to \code{student_t} with \code{df=1}), 
#'   the mean does not exist and \code{location} is the prior median.
#'   Defaults to 0, except for \code{LKJ} where there is no default, in which
#'   case how \code{location} is interpreted depends on the \code{what} 
#'   argument but always pertains to the prior \eqn{R^2} of the regression
#'   under a Beta distribution.
#' @param scale Prior scale. Default depends (see Details) but \code{Inf}
#'   implies an improper uniform prior over the positive real numbers
#' @param df Prior degrees of freedom. Defaults to 1, in which case 
#'   \code{student_t} is equivalent to \code{cauchy}.
#' @param what A character string among \code{'mode'} (the default),
#'   \code{'mean'}, \code{'median'}, or \code{'log'} indicating how the
#'   \code{location} parameter is interpred in the \code{LKJ} case. If
#'   \code{'log'}, then \code{location} is interpreted as the expected
#'   logarithm of the \eqn{R^2} under a Beta distribution. Otherwise,
#'   \code{location} is interpreted as the \code{"what"} of the \eqn{R^2}
#'   under a Beta distribution. If the number of predictors is less than
#'   or equal to two, the mode of this Beta distribution does not exist
#'   and an error will prompt the user to specify another choice for
#'   \code{"what"}.
#'   
#' @details The details depend on which prior is used
#' \subsection{Student t family}{
#'   For the prior distribution for the intercept, \code{location}, 
#'   \code{scale}, and \code{df} should be scalars. For the prior for the other
#'   coefficients they can either be vectors of length equal to the number of
#'   coefficients (not including the intercept), or they can be scalars, in 
#'   which case they will be replicated to the appropriate length.
#'   
#'   If \code{scale} is not specified it will default to 10 for the intercept
#'   and 2.5 for the other coefficients, unless the probit link function is
#'   used, in which case these defaults are scaled by a factor of 
#'   \code{dnorm(0)/dlogis(0)}, which is roughly 1.6.
#' }
#' \subsection{LKJ family}{
#'   The \code{\link{stan_lm}} and \code{\link{stan_polr}} functions allow
#'   the user to utilize a prior for the parameters called \code{LKJ}, in
#'   which case \code{location} must be a scalar on the (0,1) interval,
#'   unless \code{what = 'log'}, in which case it should be a negative
#'   scalar. The prior on all the parameters hinges on the prior beliefs about
#'   \eqn{R^2}, the proportion of variance in the outcome attributable to the
#'   predictors, which is given a \code{\link[stats]{Beta}} prior with first
#'   shape hyperparameter equal to half the number of predictors and second
#'   shape hyperparameter free. By specifying the prior mode (the default) 
#'   mean, median, or expected log of \eqn{R^2}, the second shape parameter for
#'   this Beta distribution is determined internally.
#'   
#'   For example, if \eqn{R^2 = 0.5}, then the mode, mean, and median of
#'   the \code{\link[stats]{Beta}} distribution are all the same and thus the
#'   second shape parameter is also equal to half the number of predictors.
#'   The smaller is \eqn{R^2}, the  more concentrated near zero is the prior 
#'   density for the regression coefficients. Hence, the prior on the 
#'   coefficients is regularizing and should yield a posterior distribution 
#'   with good out-of-sample predictions \emph{if} the prior is specified in a 
#'   reasonable fashion.
#' }
#' @return A named list.
#' @examples
#' \dontrun{
#' stan_glm(y ~ x1 + x2, prior = student_t(4, 0, 2.5), 
#'          prior.for.intercept = cauchy(0,10))
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
LKJ <- function(location = NULL, 
                what = c("mode", "mean", "median", "log")) {
  list(dist = "LKJ", location = location, what = what,
       df = 0, scale = 0)
}

make_eta <- function(location, what = c("mode", "mean", "median", "log"), K) {
  if (is.null(location)) stop("must specify 'location' on the (0,1) interval ",
                              "in the call to prior = LKJ(), unless 'what' is ",
                              "'log', in which case 'location' must be negative")
  stopifnot(length(location) == 1, is.numeric(location))
  stopifnot(is.numeric(K), K == as.integer(K))
  if (K == 0) stop("LKJ prior is not applicable when there are no covariates")
  what <- match.arg(what)
  half_K <- K / 2
  if (what == "mode") {
    stopifnot(location > 0, location <= 1)
    if (K <= 2) stop("mode of beta distribution does not exist when K <= 2 ",
                     "specify 'what' as 'mean', 'median', or 'log' instead")
    eta <- (half_K - 1  - location * half_K + location * 2) / location
  }
  else if (what == "mean") {
    stopifnot(location > 0, location < 1)
    eta <- (half_K - location * half_K) / location
  }
  else if (what == "median") {
    stopifnot(location > 0, location < 1)
    FUN <- function(eta) qbeta(0.5, half_K, qexp(eta)) - location
    eta <- qexp(uniroot(FUN, interval = 0:1)$root)
  }
  else { # what == "log"
    stopifnot(location < 0)
    FUN <- function(eta) digamma(half_K) - digamma(half_K + qexp(eta)) - location
    eta <- qexp(uniroot(FUN, interval = 0:1, 
                        f.lower = -location, f.upper = -.Machine$double.xmax)$root)
  }
  return(eta)
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
  