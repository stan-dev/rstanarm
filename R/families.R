# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

#' Families for rstanarm models
#' 
#' For any of the \pkg{rstanarm} modeling functions that accept a \code{family} 
#' argument, there are two ways that the family can be specified. The first way 
#' corresponds to how families are specified when calling 
#' \code{\link[stats]{glm}}, e.g., \code{family = binomial("logit")}. The second
#' method, a call to \code{rstanarm_family}, allows the user to override the
#' default values for various family-specific parameters.
#' 
#' @export
#' @param family A string naming the family object to use. Can be any of the 
#'   traditional \code{\link[stats]{family}} objects (not the "quasi" families),
#'   plus \code{\link{neg_binomial_2}} and \code{\link{t_family}}.
#' @param link A string naming the desired link function. If missing, the
#'   default link for the specified family is used.
#' @param ... Family-specific arguments related to prior distributions. The 
#'   possible arguments and their defaults are described below in Details.
#' 
#' @return A family object for use by \pkg{rstanarm} modeling functions.
#' 
#' @details
#' For all families the following two parameters can be specified in \code{...}:
#' 
#' \describe{
#'    \item{\code{min_prior_scale} (default: \code{1e-12})}{
#'      Minimum prior scale for the intercept and coefficients. The 
#'      default is nearly always fine.
#'    }
#'    \item{\code{scaled} (default: \code{TRUE})}{
#'      A logical scalar. If \code{TRUE} the prior scale is further scaled by
#'      the range of the predictor if the predictor has exactly two unique
#'      values and scaled by twice the standard deviation of the predictor if it
#'      has more than two unique values.
#'    }
#' }
#' 
#' All families except for binomial and poisson also have additional parameters
#' that can be specified in \code{...}, which are described below.
#' 
#' \subsection{Gaussian (family = "gaussian")}{
#'   \describe{
#'     \item{\code{prior_scale_for_dispersion} (default: \code{5})}{
#'       A positive scalar interpreted as the prior scale for the standard error
#'       of the regression, which is given a half-Cauchy prior truncated at 
#'       zero. This prior (and the half-Cauchy priors described for other
#'       families below) can be visualized by calling:
#'       
#'       \code{curve(2 * dcauchy(x, 0, scale = 5), from=0, to=50)}.
#'       
#'       To omit a prior (i.e., to use a flat prior) set 
#'       \code{prior_scale_for_dispersion = NULL}.
#'    }
#'   }
#' }
#' \subsection{Gamma (family = "Gamma")}{
#'   \describe{
#'     \item{\code{prior_scale_for_shape} (default: \code{5})}{
#'       A positive scalar interpreted as the prior scale for the shape
#'       parameter, which is given a half-Cauchy prior truncated at zero.
#'       To omit a prior (i.e., to use a flat prior) set 
#'       \code{prior_scale_for_shape = NULL}.
#'    }
#'   }
#' }
#' \subsection{Inverse Gaussian (family = "inverse.gaussian")}{
#'   \describe{
#'     \item{\code{prior_scale_for_shape} (default: \code{5})}{
#'       A positive scalar interpreted as the prior scale for the shape
#'       parameter, which is given a half-Cauchy prior truncated at zero.
#'       To omit a prior (i.e., to use a flat prior) set 
#'       \code{prior_scale_for_shape = NULL}.
#'    }
#'   }
#' }
#' \subsection{Negative Binomial (family = "neg_binomial_2")}{
#'   \describe{
#'     \item{\code{prior_scale_for_dispersion} (default: \code{5})}{
#'       A positive scalar interpreted as the prior scale for the overdispersion
#'       parameter, which is given a half-Cauchy prior truncated at zero.
#'       To omit a prior (i.e., to use a flat prior) set 
#'       \code{prior_scale_for_shape = NULL}.
#'    }
#'   }
#' }
#' \subsection{Student t (family = "t_family")}{
#'   \describe{
#'     \item{\code{prior_scale_for_dispersion} (default: \code{5})}{
#'       A positive scalar interpreted as the prior scale for the standard error
#'       of the regression, which is given a half-Cauchy prior truncated at
#'       zero. To omit a prior (i.e., to use a flat prior) set 
#'       \code{prior_scale_for_shape = NULL}.
#'    }
#'     \item{\code{prior_shape_for_df} (default: \code{2})}{
#'       A positive scalar interpreted as the shape parameter of a gamma prior
#'       on the degress of freedom in Student t models.
#'    }
#'    \item{\code{prior_rate_for_df} (default: \code{0.1})}{
#'      A positive scalar interpreted as the rate parameter of a gamma prior
#'      on the degress of freedom in Student t models.
#'    }
#'   }
#'   The default prior on the degrees of freedom is therefore \code{gamma(2, 
#'   0.1)}. This prior places the bulk of the prior below 30, but decays rather 
#'   slowly, as can be seen in the plot generated by 
#'   
#'   \code{curve(dgamma(x, shape = 2, rate = 0.1), from = 0, to = 100)}.
#'   
#'   To omit a prior on the degrees of freedom (i.e., to use a flat prior) set 
#'   either \code{prior_shape_for_df = NULL} or \code{prior_rate_for_df = NULL}.
#' }
#' 
#' @seealso 
#' \code{\link{neg_binomial_2}}, \code{\link{t_family}}
#' 
#' \code{\link{priors}} for specifying prior distributions for the intercept,
#' regression coefficients, covariance matrices, etc.
#' 
rstanarm_family <- function(family, link, ...) {
  family <- if (missing(link)) 
    match.fun(family)() else match.fun(family)(link = link)
  params <- list(...)
  if (length(params)) {
    # change name of any "prior_scale_for_shape" to "prior_scale_for_dispersion"
    # because that's what stan_glm.fit is expecting
    sel <- names(params) %in% "prior_scale_for_shape"
    names(params)[sel] <- "prior_scale_for_dispersion"
  }
  defaults <- default_prior_params(family)
  pars <- setdiff(names(defaults), names(params))
  if (length(pars)) {
    for (par in pars)
      params[[par]] <- defaults[[par]] 
  }
  structure(nlist(family, params), class = "rstanarm_family")
}

is.rstanarm_family <- function(x) {
  inherits(x, "rstanarm_family")
}
default_prior_params <- function(family) {
  # here family should be a family object (not just the name)
  stopifnot(is(family, "family"))
  defaults <- list(scaled = TRUE, 
                   prior_scale_for_dispersion = 5, 
                   min_prior_scale = 1e-12)
  if (is.t(family$family)) {
    defaults$prior_shape_for_df <- 2
    defaults$prior_rate_for_df <- 0.1
  } else {
    defaults$prior_shape_for_df <- 0
    defaults$prior_rate_for_df <- 0
  }
  defaults
}

#' Family function for negative binomial GLMs
#' 
#' Specifies the information required to fit a Negative Binomial GLM in a 
#' similar way to \code{\link[MASS]{negative.binomial}}. However, here the 
#' overdispersion parameter \code{theta} is not specified by the user and always
#' estimated. A call to this function can be passed to the \code{family}
#' argument of \code{\link{stan_glm}} or \code{\link{stan_glmer}} to estimate a
#' Negative Binomial model. Alternatively, the \code{\link{stan_glm.nb}} and 
#' \code{\link{stan_glmer.nb}} wrapper functions may be used, which call 
#' \code{neg_binomial_2} internally.
#' 
#' @name neg_binomial_2
#' @export
#' @param link The same as for \code{\link{poisson}}, typically a character
#'   vector of length one among \code{"log"}, \code{"identity"}, and
#'   \code{"sqrt"}.
#' @return An object of class \code{\link[stats]{family}} very similar to
#'   that of \code{\link[stats]{poisson}} but with a different family name.
#'   
#' @seealso \code{\link{rstanarm_family}} for how to set the value of the scale
#'   hyperparameter of a half-Cauchy prior on the overdispersion parameter.
#'   
#' @examples
#' fit <- stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = MASS::quine, seed = 123,
#'                 family = neg_binomial_2, QR = TRUE, 
#'                 algorithm = "fullrank") # for speed only
#' # or, equivalently, call stan_glm.nb() without specifying the family
#' 
#' # using rstanarm_family() to set prior scale for overdispersion parameter
#' nbfam <- rstanarm_family("neg_binomial_2", prior_scale_for_dispersion = 2)
#' update(fit, family = nbfam)
#'
neg_binomial_2 <- function(link = "log") {
  out <- poisson(link)
  out$family <- "neg_binomial_2"
  out$variance <- function(mu, theta) mu + mu^2 / theta
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}


#' Family function for Student t GLMs
#' 
#' @name t_family
#' @export
#' @templateVar armRef (Ch. 6)
#' @template reference-gelman-hill
#' @param link The same as for \code{\link[stats]{gaussian}}.
#' @return An object of class \code{\link[stats]{family}} very similar to
#'   that of \code{\link[stats]{gaussian}} but with a different family name.
#'   
#' @details   
#' Estimates of regression coefficients are less sensitive to outliers if the 
#' Student t distribution is used in place of the normal distribution in 
#' settings where some errors may be large. These models are sometimes referred 
#' to as \emph{robust regression models}, however that name is also used in 
#' other contexts. In the case of \pkg{rstanarm}, if a call to \code{t_family}
#' is passed to the \code{family} argument of \code{\link{stan_glm}} or
#' \code{\link{stan_glmer}} then the resulting model is a linear model where the
#' errors are believed to be distributed Student t conditional on the
#' predictors. This is equivalent to believing that the errors are Gaussian but
#' with \emph{different} standard deviations, each of which has the scaled
#' inverse chi-squared distribution.
#' 
#' These models can induce heavy-tailed posteriors which MCMC algorithms can 
#' have trouble exploring. Therefore, when fitting a model with t-distributed
#' errors, we recommend running many chains. It is also advisable to simulate
#' fake data and check that parameter values can be recovered before trusting
#' estimates when the model is fit to real data.
#' 
#' @seealso \code{\link{rstanarm_family}} for how to set the values of the
#'   hyperparameters of a gamma prior on the degrees of freedom parameter.
#'   
#' @examples 
#' SEED <- 1234
#' set.seed(SEED)
#' x <- matrix(rnorm(2000), ncol = 2)
#' alpha <- 2; beta <- c(-0.5, 0.5); df <- 4
#' y <- alpha + x %*% beta + rt(1000, df)
#' fit <- stan_glm(y ~ x, family = t_family(), seed = SEED, 
#'                 algorithm = "fullrank") # for speed only
#'                 
#' # using rstanarm_family() to set hyperparameters of gamma prior 
#' # on the degrees of freedom parameter
#' tfam <- rstanarm_family("t_family", prior_shape_for_df = 5, 
#'                         prior_rate_for_df = 0.5)
#' update(fit, family = tfam)
#' 
t_family <- function(link = "identity") {
  out <- gaussian(link)
  out$family <- "t_family"
  out$variance <- NULL
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}
