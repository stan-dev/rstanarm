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
#' @import stats
#' @import Rcpp
#' @export loo
#' @export waic
#' @export compare
#' @export launch_shinystan
#' @description An appendage to the \pkg{rstan} package that enables some of the
#'   most common applied regression models to be estimated using Markov Chain 
#'   Monte Carlo, variational approximations to the posterior distribution, or 
#'   optimization. The \pkg{rstanarm} package allows these models to be 
#'   specified using the customary R modeling syntax (e.g., like that of 
#'   \code{\link[stats]{glm}} with a \code{formula} and a \code{data.frame}).
#'   
#'   The set of models supported by \pkg{rstanarm} is large (and will continue
#'   to grow), but also limited enough so that it is possible to integrate them
#'   tightly with the \code{\link{pp_check}} function for graphical posterior
#'   predictive checks and the \code{\link{posterior_predict}} function to
#'   easily estimate the effect of specific manipulations of predictor variables
#'   or to predict the outcome in a training set. 
#'   
#'   The objects returned by the \pkg{rstanarm} modeling functions are called
#'   \code{\link[=stanreg-objects]{stanreg}} objects. In addition to all of the
#'   typical \code{\link[=stanreg-methods]{methods}} defined for fitted model
#'   objects, stanreg objects can be passed to the \code{\link[loo]{loo}}
#'   function in the \pkg{loo} package for model comparison or to the
#'   \code{\link[shinystan]{launch_shinystan}} function in the \pkg{shinystan}
#'   package in order to visualize the posterior distribution using the
#'   ShinyStan graphical user interface. See the \pkg{rstanarm} vignettes for
#'   more details about the entire process.
#'
#' @section Estimation algorithms: 
#' The modeling functions in the \pkg{rstanarm} package take an \code{algorithm}
#' argument that can be one of the following:
#' \describe{
#'  \item{\strong{Sampling} (\code{algorithm="sampling"})}{
#'  Uses Markov Chain Monte Carlo (MCMC) --- in particular, Hamiltonian Monte 
#'  Carlo (HMC) with a tuned but diagonal mass matrix --- to draw from the 
#'  posterior distribution of the parameters. See \code{\link[rstan]{sampling}} 
#'  for more details. This is the slowest but most reliable of the available
#'  estimation algorithms and it is \strong{the default and recommended
#'  algorithm for statistical inference.}
#'  }
#'  \item{\strong{Mean-field} (\code{algorithm="meanfield"})}{
#'  Uses mean-field variational inference to draw from an approximation to the
#'  posterior distribution. In particular, this algorithm finds the set of
#'  independent normal distributions in the unconstrained space that --- when
#'  transformed into the constrained space --- most closely approximate the
#'  posterior distribution. Then it draws repeatedly from these independent
#'  normal distributions and transforms them into the constrained space. The
#'  entire process is much faster than HMC and yields independent draws but
#'  \strong{is not recommended for final statistical inference}. It can be
#'  useful to narrow the set of candidate models in large problems, particularly
#'  when specifying \code{QR=TRUE} in \code{\link{stan_glm}},
#'  \code{\link{stan_glmer}}, and \code{\link{stan_gamm4}}, but is \strong{only
#'  an approximation to the posterior distribution}.
#'  }
#'  \item{\strong{Full-rank} (\code{algorithm="fullrank"})}{
#'  Uses full-rank variational inference to draw from an approximation to the 
#'  posterior distribution by finding the multivariate normal distribution in 
#'  the unconstrained space that --- when transformed into the constrained space
#'  --- most closely approximates the posterior distribution. Then it draws 
#'  repeatedly from this multivariate normal distribution and transforms the 
#'  draws into the constrained space. This process is slower than meanfield 
#'  variational inference but is faster than HMC. Although still an 
#'  approximation to the posterior distribution and thus \strong{not recommended
#'  for final statistical inference}, the approximation is more realistic than 
#'  that of mean-field variational inference because the parameters are not 
#'  assumed to be independent in the unconstrained space. Nevertheless, fullrank
#'  variational inference is a more difficult optimization problem and the 
#'  algorithm is more prone to non-convergence or convergence to a local 
#'  optimum.
#'  }
#'  \item{\strong{Optimizing} (\code{algorithm="optimizing"})}{
#'  Finds the posterior mode using a C++ implementation of the LBGFS algorithm.
#'  See \code{\link[rstan]{optimizing}} for more details. If there is no prior 
#'  information, then this is equivalent to maximum likelihood, in which case 
#'  there is no great reason to use the functions in the \pkg{rstanarm} package 
#'  over the emulated functions in other packages. However, if priors are 
#'  specified, then the estimates are penalized maximum likelihood estimates, 
#'  which may have some redeeming value. Currently, optimization is only 
#'  supported for \code{\link{stan_glm}}.
#'  }
#' }
#' 
#' 
#' @section Modeling functions: 
#' The model estimating functions are described in greater detail in their
#' individual help pages and vignettes. Here we provide a very brief
#' overview:
#' 
#' \describe{
#'  \item{\code{\link{stan_lm}}, \code{stan_aov}}{
#'   Similar to \code{\link[stats]{lm}} or \code{\link[stats]{aov}} but with 
#'   novel regularizing priors on the model parameters that are driven by prior 
#'   beliefs about \eqn{R^2}, the proportion of variance in the outcome 
#'   attributable to the predictors in a linear model.
#'  }
#'  \item{\code{\link{stan_glm}}, \code{stan_glm.nb}}{
#'   Similar to \code{\link[stats]{glm}} but with Gaussian, Student t, Cauchy 
#'   or hierarhical shrinkage prior distributions for the coefficients and,
#'   if applicable, a half-Cauchy prior for any nuisance parameter in a 
#'   Generalized Linear Model (GLM) that is characterized by a 
#'   \code{\link[stats]{family}} object. It is also possible to estimate a 
#'   negative bionomial model in a similar way to the \code{\link[MASS]{glm.nb}} 
#'   function in the \pkg{MASS} package.
#'  }
#'  \item{\code{\link{stan_glmer}}, \code{stan_glmer.nb}, \code{stan_lmer}}{
#'   Similar to the \code{\link[lme4]{glmer}}, \code{\link[lme4]{glmer.nb}} and 
#'   \code{\link[lme4]{lmer}} functions in the \pkg{lme4} package in that GLMs 
#'   are augmented to have group-specific terms that deviate from the common 
#'   coefficients according to a mean-zero multivariate normal distribution with
#'   a highly-structured but unknown covariance matrix (for which \pkg{rstanarm}
#'   introduces an innovative prior distribution). MCMC provides more
#'   appropriate estimates of uncertainty for models that consist of a mix of
#'   common and group-specific parameters.
#'  }
#'  \item{\code{\link{stan_gamm4}}}{
#'   Similar to \code{\link[gamm4]{gamm4}} in the \pkg{gamm4} package, which 
#'   augments a GLM (possibly with group-specific terms) with nonlinear smooth 
#'   functions of the predictors to form a Generalized Additive Mixed Model 
#'   (GAMM). Rather than calling \code{\link[lme4]{glmer}} like 
#'   \code{\link[gamm4]{gamm4}} does, \code{\link{stan_gamm4}} essentially calls
#'   \code{\link{stan_glmer}}, which avoids the optimization issues that often 
#'   crop up with GAMMs and provides better estimates for the uncertainty of the
#'   parameter estimates.
#'  }
#'  \item{\code{\link{stan_polr}}}{
#'   Similar to \code{\link[MASS]{polr}} in the \pkg{MASS} package in that it
#'   models an ordinal response but also implies a prior distribution on the 
#'   unknown cutpoints. Can also be used to model binary outcomes, possibly
#'   while estimating an unknown exponent governing the probability of success.
#'  }
#' }
#' 
#' @section Prior distributions:
#' See \code{\link{priors}} for an overview of the various choices the user can 
#' make for prior distributions. The package vignettes also provide 
#' examples of using many of the available priors as well as more detailed 
#' descriptions of some of the novel priors used by \pkg{rstanarm}.
#'  
#' @seealso \code{\link{stanreg-objects}} and \code{\link{stanreg-methods}} for 
#'   details on the fitted model objects returned by the modeling functions.
#'   
#'   \code{\link{plots}} for the various plots that can be used
#'   to explore and check fitted models.
#'   
#'   \url{http://mc-stan.org/} for more information on the Stan C++ package used
#'   by \pkg{rstanarm} for model fitting.
#'   
#'   \url{https://github.com/stan-dev/rstanarm/issues/} to submit a bug
#'   report or feature request.
#'   
#'   \url{https://groups.google.com/forum/#!forum/stan-users/} to ask a question 
#'   about \pkg{rstanarm} on the Stan-users forum.
#'   
#' @templateVar armRef \url{http://stat.columbia.edu/~gelman/arm/}
#' @templateVar bdaRef \url{http://stat.columbia.edu/~gelman/book/}
#' @template reference-lme4
#' @template reference-bda
#' @template reference-gelman-hill
#' @template reference-stan-manual
#' @template reference-loo
#'   
NULL
