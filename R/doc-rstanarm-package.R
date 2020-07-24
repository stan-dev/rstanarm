# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

#' Applied Regression Modeling via RStan
#'
#' @docType package
#' @name rstanarm-package
#' @aliases rstanarm
#' @useDynLib rstanarm, .registration = TRUE
#'
#' @import methods
#' @importFrom rstan optimizing sampling vb constrain_pars extract
#'   extract_sparse_parts get_posterior_mean stanc
#' @importFrom utils capture.output
#' @import stats
#' @import Rcpp
#' @import bayesplot
#' @import shinystan
#' @import rstantools
#' @export log_lik posterior_linpred posterior_epred posterior_predict posterior_interval
#' @export predictive_interval predictive_error prior_summary bayes_R2
#' @export loo_linpred loo_predict loo_predictive_interval loo_R2
#' @export loo waic kfold loo_compare
#' @export launch_shinystan
#'
#' @description
#' \if{html}{
#'    \figure{stanlogo.png}{options: width="50px" alt="http://mc-stan.org/about/logo/"}
#'    \emph{Stan Development Team}
#' }
#'
#' The \pkg{rstanarm} package is an appendage to the \pkg{rstan} package that
#' enables many of the most common applied regression models to be estimated
#' using Markov Chain Monte Carlo, variational approximations to the posterior
#' distribution, or optimization. The \pkg{rstanarm} package allows these models
#' to be specified using the customary R modeling syntax (e.g., like that of
#' \code{\link[stats]{glm}} with a \code{formula} and a \code{data.frame}).
#'
#' The sections below provide an overview of the modeling functions and
#' estimation algorithms used by \pkg{rstanarm}.
#'
#' @details
#' The set of models supported by \pkg{rstanarm} is large (and will continue to
#' grow), but also limited enough so that it is possible to integrate them
#' tightly with the \code{\link{pp_check}} function for graphical posterior
#' predictive checks with \pkg{\link[bayesplot:bayesplot-package]{bayesplot}} and the
#' \code{\link{posterior_predict}} function to easily estimate the effect of
#' specific manipulations of predictor variables or to predict the outcome in a
#' training set.
#'
#' The objects returned by the \pkg{rstanarm} modeling functions are called
#' \code{\link[=stanreg-objects]{stanreg}} objects. In addition to all of the
#' typical \code{\link[=stanreg-methods]{methods}} defined for fitted model
#' objects, stanreg objects can be passed to the \code{\link[loo]{loo}} function
#' in the \pkg{loo} package for model comparison or to the
#' \code{\link[shinystan]{launch_shinystan}} function in the \pkg{shinystan}
#' package in order to visualize the posterior distribution using the ShinyStan
#' graphical user interface. See the \pkg{rstanarm} vignettes for more details
#' about the entire process.
#'
#' @inheritSection available-models Modeling functions
#' @inheritSection available-algorithms Estimation algorithms
#'
#' @section Prior distributions:
#' See \link[=priors]{priors help page} and the vignette
#' \href{http://mc-stan.org/rstanarm/articles/priors.html}{\emph{Prior Distributions for rstanarm Models}}
#' for an overview of the various choices the user can make for prior
#' distributions. The package vignettes for the modeling functions also provide
#' examples of using many of the available priors as well as more detailed
#' descriptions of some of the novel priors used by \pkg{rstanarm}.
#'
#' @seealso
#' \itemize{
#'   \item \url{http://mc-stan.org/} for more information on the Stan C++
#'   package used by \pkg{rstanarm} for model fitting.
#'   \item \url{https://github.com/stan-dev/rstanarm/issues/} to submit a bug
#'   report or feature request.
#'   \item \url{http://discourse.mc-stan.org} to ask a
#'   question about \pkg{rstanarm} on the Stan-users forum.
#' }
#'
#' @templateVar armRef \url{http://stat.columbia.edu/~gelman/arm/}
#' @templateVar bdaRef \url{http://stat.columbia.edu/~gelman/book/}
#' @template reference-lme4
#' @template reference-bda
#' @template reference-gelman-hill
#' @template reference-stan-manual
#' @template reference-loo
#' @template reference-bayesvis
#' @template reference-muth
#'
NULL
