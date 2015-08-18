#' Prior distributions
#' 
#' These functions are used to specify the prior-related arguments of various
#' model-fitting functions in the \pkg{rstanarm} package.
#' 
#' @export 
#' @name priors
#' @param location Prior location. For \code{normal} and 
#'   \code{student_t} (provided that \code{df > 1}) this is the prior mean. For 
#'   \code{cauchy} (which is equivalent to \code{student_t} with \code{df=1}), 
#'   the mean does not exist and \code{location} is the prior median.
#'   Defaults to 0, except for \code{R2} where there is no default, in which
#'   case how \code{location} is interpreted depends on the \code{what} 
#'   argument but always pertains to the prior location of the \eqn{R^2} 
#'   under a Beta distribution. See the Details section.
#' @param scale Prior scale. Default depends on the family (see Details)
#' @param df,df1,df2 Prior degrees of freedom. Defaults to 1, in which case 
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
#' @details The details depend on the family of the prior being used:
#' \subsection{Student t family}{
#'   For the prior distribution for the intercept, \code{location}, 
#'   \code{scale}, and \code{df} should be scalars. For the prior for the other
#'   coefficients they can either be vectors of length equal to the number of
#'   coefficients (not including the intercept), or they can be scalars, in 
#'   which case they will be recycled to the appropriate length. As the 
#'   degrees of freedom approaches infinity, the Student t distribution 
#'   approaches the normal distribution and if the degrees of freedom are one,
#'   then the Student t distribution is the Cauchy distribution.
#'   
#'   If \code{scale} is not specified it will default to 10 for the intercept
#'   and 2.5 for the other coefficients, unless the probit link function is
#'   used, in which case these defaults are scaled by a factor of 
#'   \code{dnorm(0)/dlogis(0)}, which is roughly 1.6.
#' }
#' \subsection{Hierarchical shrinkage family}{
#'   The horseshoe prior is a normal with a mean of zero and a standard 
#'   deviation that is distributed half Student t with some 
#'   (traditionally 1 but by default 3) degrees of freedom, scaled by a
#'   half Cauchy parameter.
#'   
#'   The horseshoe plus prior is a normal with a mean of zero and a standard
#'   deviation that is distributed as the product of two independent half 
#'   Student t parameters with some (traditionally 1 but by default 3) 
#'   degrees of freedom that are each scaled by the same square root of a
#'   half Cauchy parameter.
#'   
#'   These hierarchical shrinkage priors have very tall modes and very fat
#'   tails. Consequently, the tend to produce posterior distributions that
#'   are very concentrated near zero, unless the predictor has a strong
#'   influence on the outcome, in which case the prior has little influence.
#'   Traditionally, horseshoe priors set the degrees of freedom equal to 1,
#'   but that can make it difficult for Stan to sample without encountering
#'   some divergent transitions. Thus, the default degrees of freedom for
#'   \code{horseshoe} and \code{horseshoe_plus} is 3, which makes the 
#'   variance of the Student t distribution finite and produces some shrinkage
#'   on the coefficients, even for strong predictors.
#' }
#' \subsection{Covariance matrices}{
#'   Covariance matrices are decomposed into correlation matrices and 
#'   variances, and the variances are in turn decomposed into a scale parameter
#'   and a simplex vector. This prior is represented by the \code{decov} 
#'   function.
#'   
#'   The prior for a correlation matrix is called LKJ whose density is 
#'   proportional to the determinant of the correlation matrix raised to the 
#'   power of \eqn{shape - 1}, which depends solely on the positive shape 
#'   parameter. If \code{shape = 1} (the default), then this prior is jointly 
#'   uniform over all correlation matrices of that size. If \code{shape > 1},
#'   then the identity matrix is the mode (it is always the mean) and in the 
#'   unlikely case that \code{shape < 1}, the identity matrix is the trough.
#'   
#'   The trace of a covariance matrix is equal to the sum of the variances and
#'   we set the trace equal to the product of the size of the covariance matrix
#'   and the \emph{square} of a positive \code{scale} parameter. The particular
#'   variances are set equal to the product of a simplex vector --- which is
#'   non-negative and sums to \eqn{1} --- and the scalar trace. In other words,
#'   each element of the simplex vector represents the proportion of the trace
#'   attributable to the corresponding variable (the stick-breaking metaphor).
#'   
#'   A symmetric Dirichlet prior is used for a simplex vector, which has a 
#'   single (positive) \code{concentration} parameter, which defaults to
#'   \eqn{1} and implies that the prior is jointly uniform over the space of
#'   simplex vectors of that size. If \code{concentration > 1}, then the prior
#'   mode corresponds to all variables having the same (proportion of total)
#'   variance, which can be used to ensure the the posterior variances are not
#'   zero. As the \code{concentration} parameter approaches infinity, this
#'   mode becomes more pronounced. In the unlikely case that 
#'   \code{concentration < 1}, the variances are more polarized.
#'   
#'   If all the variables were multiplied by a number, the trace of their
#'   covariance matrix would increase by that number squared. Thus, it is
#'   reasonable to use a scale-invariant prior distribution, and in this case 
#'   we utilize a gamma distribution, whose \code{shape_gamma} and \code{scale} 
#'   parameters are both \eqn{1} by default, implying a unit-exponential 
#'   distribution. We scale up by the square root of the number of variables to 
#'   make the default value of the \code{scale} parameter more widely 
#'   applicable. Set the \code{shape_gamma} hyperparameter to some value
#'   greater than one to ensure that the posterior trace is not zero.
#'   
#'   If \code{shape}, \code{concentration}, \code{shape_gamma} and / or 
#'   \code{scale} are positive scalars, then they are recycled to the 
#'   appropriate length. Otherwise, each can be a positive vector of the 
#'   appropriate length, but the appropriate length depends on the number of 
#'   covariance matrices in the model and their sizes. A one-by-one covariance 
#'   matrix is just a variance and thus does not have \code{shape} or 
#'   \code{concentration} parameters, but does have \code{shape_gamma} and 
#'   \code{scale} parameter for the the prior standard deviation of that
#'   variable.
#' }
#' \subsection{R2 family}{
#'   The \code{\link{stan_lm}}, \code{\link{stan_aov}} and 
#'   \code{\link{stan_polr}} functions allow the user to utilize a function 
#'   called \code{R2} to convey prior information about all the parameters. 
#'   This prior hinges on prior beliefs about the location of \eqn{R^2}, the 
#'   proportion of variance in the outcome attributable to the predictors, 
#'   which has a \code{\link[stats]{Beta}} prior with first shape 
#'   hyperparameter equal to half the number of predictors and second shape 
#'   hyperparameter free. By specifying the prior mode (the default) mean, 
#'   median, or expected log of \eqn{R^2}, the second shape parameter for this 
#'   Beta distribution is determined internally. If \code{what = 'log'}, 
#'   location should be a negative scalar; otherwise it should be a scalar on 
#'   the \eqn{(0,1)} interval.
#'   
#'   For example, if \eqn{R^2 = 0.5}, then the mode, mean, and median of
#'   the \code{\link[stats]{Beta}} distribution are all the same and thus the
#'   second shape parameter is also equal to half the number of predictors.
#'   The second shape parameter of the \code{\link[stats]{Beta}} distribution
#'   is actually the same as the shape parameter in the LKJ prior for a
#'   correlation matrix described in the previous subsection. Thus, the smaller 
#'   is \eqn{R^2}, the larger is the shape parameter, the smaller are the
#'   prior correlations among the outcome and predictor variables, and the more
#'   concentrated near zero is the prior density for the regression 
#'   coefficients. Hence, the prior on the coefficients is regularizing and 
#'   should yield a posterior distribution with good out-of-sample predictions 
#'   \emph{if} the prior location of \eqn{R^2} is specified in a reasonable 
#'   fashion.
#' }
#' @return A named list.
#' @examples
#' stan_glm(mpg ~ ., data = mtcars, prior_PD = TRUE, chains = 1,
#'          prior = student_t(4, 0, 2.5), prior.for.intercept = cauchy(0,10), 
#'          prior.options = prior_options(prior.scale.for.dispersion = 2))
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
horseshoe <- function(df = 3) {
  validate_parameter_value(df)
  nlist(dist = "horseshoe", df, location = 0, scale = 1)
}

#' @rdname priors
#' @export
horseshoe_plus <- function(df1 = 3, df2 = 3) {
  validate_parameter_value(df1)
  validate_parameter_value(df2)
  # scale gets used as a second df hyperparameter
  nlist(dist = "horseshoe_plus", df = df1, location = 0, scale = df2)
}

#' @rdname priors
#' @export
cauchy <- function(location = 0, scale = NULL) {
  student_t(df = 1, location = location, scale = scale)
}

#' @rdname priors
#' @param shape Shape parameter for an LKJ prior on the correlation matrix 
#'  in the \code{decov} prior.
#' @param concentration Concentration parameter for the symmetric Dirichlet 
#'  distribution in the \code{decov} prior.
#' @param gamma_shape Shape parameter for a gamma prior on the scale 
#'   parameter in the \code{dcov} prior.

#' @export
decov <- function(shape = 1, concentration = 1, gamma_shape = 1, scale = 1) {
  if (any(shape <= 0)) stop("'shape' parameter must be positive")
  if (any(concentration <= 0)) stop("'concentration' parameter must be positive")
  if (any(gamma_shape <= 0)) stop("'gamma_shape' parameter must be positive")
  if (any(scale <= 0)) stop("'scale' parameter must be positive")
  nlist(shape, concentration, gamma_shape, scale)
}

#' @rdname priors
#' @export
R2 <- function(location = NULL, 
                what = c("mode", "mean", "median", "log")) {
  list(dist = "R2", location = location, what = what,
       df = 0, scale = 0)
}

make_eta <- function(location, what = c("mode", "mean", "median", "log"), K) {
  if (is.null(location)) stop("must specify 'location' on the (0,1) interval ",
                              "in the call to prior = R2(), unless 'what' is ",
                              "'log', in which case 'location' must be negative")
  stopifnot(length(location) == 1, is.numeric(location))
  stopifnot(is.numeric(K), K == as.integer(K))
  if (K == 0) stop("R2 prior is not applicable when there are no covariates")
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
#' @param scaled A logical scalar, defaulting to \code{TRUE}. If \code{TRUE} the
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
  