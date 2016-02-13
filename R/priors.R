# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Prior distributions and options
#' 
#' These functions are used to specify the prior-related arguments of various
#' modeling functions in the \pkg{rstanarm} package.
#' 
#' @export 
#' @name priors
#' @param location Prior location. For \code{normal} and \code{student_t} 
#'   (provided that \code{df > 1}) this is the prior mean. For \code{cauchy} 
#'   (which is equivalent to \code{student_t} with \code{df=1}), the mean does 
#'   not exist and \code{location} is the prior median. The default value is 
#'   \eqn{0}, except for \code{R2} which has no default value for
#'   \code{location}. For \code{R2}, \code{location} pertains to the prior
#'   location of the \eqn{R^2} under a Beta distribution, but the interpretation
#'   of the \code{location} parameter depends on the specified value of the
#'   \code{what} argument (see Details).
#' @param scale Prior scale. The default depends on the family (see Details).
#' @param df,df1,df2 Prior degrees of freedom. The default is \eqn{1} for 
#'   \code{student_t}, in which case it is equivalent to \code{cauchy}. For the
#'   hierarchical shrinkage priors (\code{hs} and \code{hs_plus}) the degrees of
#'   freedom parameter(s) default to \eqn{3}.
#' @param what A character string among \code{'mode'} (the default),
#'   \code{'mean'}, \code{'median'}, or \code{'log'} indicating how the
#'   \code{location} parameter is interpreted in the \code{LKJ} case. If
#'   \code{'log'}, then \code{location} is interpreted as the expected
#'   logarithm of the \eqn{R^2} under a Beta distribution. Otherwise,
#'   \code{location} is interpreted as the \code{what} of the \eqn{R^2}
#'   under a Beta distribution. If the number of predictors is less than
#'   or equal to two, the mode of this Beta distribution does not exist
#'   and an error will prompt the user to specify another choice for
#'   \code{what}.
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
#'   The hierarchical shrinkage priors are normal with a mean of zero and a 
#'   standard deviation that is also a random variable. The traditional 
#'   hierarchical shrinkage prior utilizes a standard deviation that is 
#'   distributed half Cauchy with a median of zero and a scale parameter that is
#'   also half Cauchy. This is called the "horseshoe prior". The hierarchical 
#'   shrinkage (\code{hs}) prior in the \pkg{rstanarm} package instead utilizes 
#'   a half Student t distribution for the standard deviation (with 3 degrees of
#'   freedom by default), scaled by a half Cauchy parameter, as described by
#'   Piironen and Vehtari (2015). It is possible to change the \code{df}
#'   argument, the prior degrees of freedom, to obtain less or more shrinkage.
#'   
#'   The hierarhical shrinkpage plus (\code{hs_plus}) prior is a normal with a 
#'   mean of zero and a standard deviation that is distributed as the product of
#'   two independent half Student t parameters (both with 3 degrees of freedom
#'   (\code{df1}, \code{df2}) by default) that are each scaled by the same
#'   square root of a half Cauchy parameter.
#'   
#'   These hierarchical shrinkage priors have very tall modes and very fat 
#'   tails. Consequently, they tend to produce posterior distributions that are
#'   very concentrated near zero, unless the predictor has a strong influence on
#'   the outcome, in which case the prior has little influence. Hierarchical
#'   shrinkage priors often require you to increase the 
#'   \code{\link{adapt_delta}} tuning parameter in order to diminish the number
#'   of divergent transitions. For more details on tuning parameters and
#'   divergent transitions see the Troubleshooting section of the 
#'   \emph{How to Use the rstanarm Package} vignette.
#' }
#' \subsection{Dirichlet family}{
#'   The Dirichlet distribution is a multivariate generalization of the beta
#'   distribution. It is perhaps the easiest prior distribution to specify
#'   because the concentration parameters can be interpreted as prior counts
#'   (although they need not be integers) of a multinomial random variable.
#'   
#'   The Dirichlet distribution is used in \code{\link{stan_polr}} for an 
#'   implicit prior on the cutpoints in an ordinal regression model. More
#'   specifically, the Dirichlet prior pertains to the prior probability of
#'   observing each category of the ordinal outcome when the predictors are at
#'   their sample means. Given these prior probabilities, it is straightforward
#'   to add them to form cumulative probabilities and then use an inverse CDF
#'   transformation of the cumulative probabilities to define the cutpoints.
#'   
#'   If a scalar is passed to the \code{concentration} argument of the 
#'   \code{dirichlet} function, then it is replicated to the appropriate length 
#'   and the Dirichlet distribution is symmetric. If \code{concentration} is a
#'   vector and all elements are \eqn{1}, then the Dirichlet distribution is
#'   jointly uniform. If all concentration parameters are equal but greater than
#'   \eqn{1} then the prior mode is that the categories are equiprobable, and
#'   the larger the value of the identical concentration parameters, the more
#'   sharply peaked the distribution is at the mode. The elements in 
#'   \code{concentration} can also be given different values to represent that 
#'   not all outcome categories are a priori equiprobable.
#' }
#' \subsection{Covariance matrices}{
#'   Covariance matrices are decomposed into correlation matrices and 
#'   variances. The variances are in turn decomposed into the product of a
#'   simplex vector and the trace of the matrix. Finally, the trace is the
#'   product of the order of the matrix and the square of a scale parameter.
#'   This prior on a covariance matrix is represented by the \code{decov} 
#'   function.
#'   
#'   The prior for a correlation matrix is called LKJ whose density is 
#'   proportional to the determinant of the correlation matrix raised to the 
#'   power of a positive regularization parameter minus one. If
#'   \code{regularization = 1} (the default), then this prior is jointly 
#'   uniform over all correlation matrices of that size. If 
#'   \code{regularization > 1}, then the identity matrix is the mode and in the
#'   unlikely case that \code{regularization < 1}, the identity matrix is the
#'   trough.
#'   
#'   The trace of a covariance matrix is equal to the sum of the variances. We
#'   set the trace equal to the product of the order of the covariance matrix
#'   and the \emph{square} of a positive scale parameter. The particular
#'   variances are set equal to the product of a simplex vector --- which is
#'   non-negative and sums to \eqn{1} --- and the scalar trace. In other words,
#'   each element of the simplex vector represents the proportion of the trace
#'   attributable to the corresponding variable.
#'   
#'   A symmetric Dirichlet prior is used for the simplex vector, which has a 
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
#'   reasonable to use a scale-invariant prior distribution for the positive
#'   scale parameter, and in this case we utilize a Gamma distribution, whose
#'   \code{shape} and \code{scale} are both \eqn{1} by default, implying a
#'   unit-exponential distribution. Set the \code{shape} hyperparameter to some
#'   value greater than \eqn{1} to ensure that the posterior trace is not zero.
#'   
#'   If \code{regularization}, \code{concentration}, \code{shape} and / or 
#'   \code{scale} are positive scalars, then they are recycled to the 
#'   appropriate length. Otherwise, each can be a positive vector of the 
#'   appropriate length, but the appropriate length depends on the number of 
#'   covariance matrices in the model and their sizes. A one-by-one covariance 
#'   matrix is just a variance and thus does not have \code{regularization} or 
#'   \code{concentration} parameters, but does have \code{shape} and 
#'   \code{scale} parameters for the prior standard deviation of that 
#'   variable.
#' }
#' \subsection{R2 family}{
#'   The \code{\link{stan_lm}}, \code{\link{stan_aov}}, and 
#'   \code{\link{stan_polr}} functions allow the user to utilize a function 
#'   called \code{R2} to convey prior information about all the parameters. 
#'   This prior hinges on prior beliefs about the location of \eqn{R^2}, the 
#'   proportion of variance in the outcome attributable to the predictors, 
#'   which has a \code{\link[stats]{Beta}} prior with first shape 
#'   hyperparameter equal to half the number of predictors and second shape 
#'   hyperparameter free. By specifying the prior mode (the default), mean, 
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
#' @return A named list to be used internally by the \pkg{rstanarm} model
#'   fitting functions.
#' @seealso The various vignettes for the \pkg{rstanarm} package also discuss 
#'   and demonstrate the use of some of the supported prior distributions.
#' 
#' @templateVar bdaRef \url{http://stat.columbia.edu/~gelman/book/}
#' @template reference-bda
#' @template reference-piironen-vehtari
#' @template reference-stan-manual
#' 
#' @examples
#' fmla <- mpg ~ wt + qsec + drat + am
#' 
#' # Draw from prior predictive distribution (by setting prior_PD = TRUE)
#' prior_pred_fit <- stan_glm(fmla, data = mtcars, chains = 1, prior_PD = TRUE,
#'                            prior = student_t(df = 4, 0, 2.5), 
#'                            prior_intercept = cauchy(0,10), 
#'                            prior_ops = prior_options(prior_scale_for_dispersion = 2))
#' 
#' \dontrun{
#' # Can assign priors to names
#' N05 <- normal(0, 5)
#' fit <- stan_glm(fmla, data = mtcars, prior = N05, prior_intercept = N05)
#' }
#' 
#' # Visually compare normal, student_t, and cauchy
#' library(ggplot2)
#' compare_priors <- function(scale = 1, df_t = 2, xlim = c(-10, 10)) {
#'   dt_loc_scale <- function(x, df, location, scale) { 
#'     # t distribution with location & scale parameters
#'     1 / scale * dt((x - location) / scale, df)  
#'   }
#'   ggplot(data.frame(x = xlim), aes(x)) + 
#'     stat_function(fun = dnorm, 
#'                   args = list(mean = 0, sd = scale), 
#'                   color = "purple", size = .75) +
#'     stat_function(fun = dt_loc_scale, 
#'                   args = list(df = df_t, location = 0, scale = scale), 
#'                   color = "orange", size = .75) +
#'     stat_function(fun = dcauchy, 
#'                   args = list(location = 0, scale = scale), 
#'                   color = "skyblue", size = .75, linetype = 2) + 
#'     ggtitle("normal (purple) vs student_t (orange) vs cauchy (blue)")
#' }
#' # Cauchy has fattest tails, then student_t, then normal
#' compare_priors()
#' 
#' # The student_t with df = 1 is the same as the cauchy
#' compare_priors(df_t = 1) 
#' 
#' # Even a scale of 5 is somewhat large. It gives plausibility to rather 
#' # extreme values
#' compare_priors(scale = 5, xlim = c(-20,20)) 
#' 
#' # If you use a prior like normal(0, 1000) to be "non-informative" you are 
#' # actually saying that a coefficient value of e.g. -500 is quite plausible
#' compare_priors(scale = 1000, xlim = c(-1000,1000))
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
hs <- function(df = 3) {
  validate_parameter_value(df)
  nlist(dist = "hs", df, location = 0, scale = 1)
}

#' @rdname priors
#' @export
hs_plus <- function(df1 = 3, df2 = 3) {
  validate_parameter_value(df1)
  validate_parameter_value(df2)
  # scale gets used as a second df hyperparameter
  nlist(dist = "hs_plus", df = df1, location = 0, scale = df2)
}

#' @rdname priors
#' @export
#' @param regularization Exponent for an LKJ prior on the correlation matrix in
#'   the \code{decov} prior. The default is \eqn{1}, implying a joint uniform
#'   prior.
#' @param concentration Concentration parameter for a symmetric Dirichlet 
#'   distribution. The defaults is \eqn{1}, implying a joint uniform prior.
#' @param shape Shape parameter for a gamma prior on the scale parameter in the
#'   \code{decov} prior. If \code{shape} and \code{scale} are both \eqn{1} (the
#'   default) then the gamma prior simplifies to the unit-exponential
#'   distribution.
decov <- function(regularization = 1, concentration = 1, 
                  shape = 1, scale = 1) {
  validate_parameter_value(regularization)
  validate_parameter_value(concentration)
  validate_parameter_value(shape)
  validate_parameter_value(scale)
  nlist(dist = "decov", regularization, concentration, shape, scale)
}

#' @rdname priors
#' @export
dirichlet <- function(concentration = 1) {
  validate_parameter_value(concentration)
  nlist(dist = "dirichlet", concentration)
}

#' @rdname priors
#' @export
R2 <- function(location = NULL, what = c("mode", "mean", "median", "log")) {
  list(dist = "R2", location = location, what = what, df = 0, scale = 0)
}

#' @rdname priors
#' @export 
#' @param prior_scale_for_dispersion Prior scale for the standard error of the 
#'   regression in Gaussian models, which is given a half-Cauchy prior truncated
#'   at zero.
#' @param min_prior_scale Minimum prior scale for the intercept and 
#'   coefficients.
#' @param scaled A logical scalar, defaulting to \code{TRUE}. If \code{TRUE} the
#'   \code{prior_scale} is further scaled by the range of the predictor if the 
#'   predictor has exactly two unique values and scaled by twice the standard
#'   deviation of the predictor if it has more than two unique values.
#'
prior_options <- function(prior_scale_for_dispersion = 5, 
                          min_prior_scale = 1e-12, 
                          scaled = TRUE) {
  validate_parameter_value(prior_scale_for_dispersion)
  validate_parameter_value(min_prior_scale)
  nlist(scaled, min_prior_scale, prior_scale_for_dispersion)
}


make_eta <- function(location, what = c("mode", "mean", "median", "log"), K) {
  if (is.null(location)) 
    stop("For the R2 prior, 'location' must be in the (0,1) interval unless ",
         "'what' is 'log'. If 'what' is 'log' then 'location' must be negative.",
         call. = FALSE)
  stopifnot(length(location) == 1, is.numeric(location))
  stopifnot(is.numeric(K), K == as.integer(K))
  if (K == 0) 
    stop("R2 prior is not applicable when there are no covariates.", 
         call. = FALSE)
  what <- match.arg(what)
  half_K <- K / 2
  if (what == "mode") {
    stopifnot(location > 0, location <= 1)
    if (K <= 2)
      stop(paste("R2 prior error.", 
                 "The mode of the beta distribution does not exist",
                 "with fewer than three predictors.", 
                 "Specify 'what' as 'mean', 'median', or 'log' instead."),
           call. = FALSE)
    eta <- (half_K - 1  - location * half_K + location * 2) / location
  } else if (what == "mean") {
    stopifnot(location > 0, location < 1)
    eta <- (half_K - location * half_K) / location
  } else if (what == "median") {
    stopifnot(location > 0, location < 1)
    FUN <- function(eta) qbeta(0.5, half_K, qexp(eta)) - location
    eta <- qexp(uniroot(FUN, interval = 0:1)$root)
  } else { # what == "log"
    stopifnot(location < 0)
    FUN <- function(eta) digamma(half_K) - digamma(half_K + qexp(eta)) - location
    eta <- qexp(uniroot(FUN, interval = 0:1, 
                        f.lower = -location, 
                        f.upper = -.Machine$double.xmax)$root)
  }
  
  return(eta)
}
