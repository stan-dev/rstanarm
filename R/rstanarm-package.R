#' Applied Regression Modeling via RStan
#'
#' @docType package
#' @name rstanarm-package
#' @aliases rstanarm
#' @useDynLib rstanarm, .registration = TRUE 
#' 
#' @import rstan
#' @import stats
#' @import Rcpp
#' 
#' @export loo
#' @export waic
#' @export launch_shinystan
#' 
#' @description An appendage to the \pkg{rstan} package that enables the most
#'   common applied regression models to be estimated using Markov Chain Monte
#'   Carlo (or optimization) but still specified using customary R modeling
#'   syntax (e.g. like that of \code{\link[stats]{glm}} with a
#'   \code{\link[stats]{formula}} and a \code{\link{data.frame}}).
#'   
#'   The set of models supported by the \pkg{rstanarm} package is extensive (and
#'   will continue to grow), but also limited enough so that it is possible to 
#'   integrate them tightly with the \code{\link{posterior_predict}} function to
#'   estimate the effect of specific manipulations of predictor variables or to 
#'   predict the outcome in a training set. Also, any 
#'   \code{\link[=stanreg-objects]{stanreg}} objects can be passed to the 
#'   \code{\link[loo]{loo}} function in the \pkg{loo} package for model 
#'   comparison or to the \code{\link[shinystan]{launch_shinystan}} function in 
#'   the \pkg{shinystan} package in order to better visualize the posterior 
#'   distribution using the ShinyStan graphical user interface. See the
#'   \pkg{rstanarm} vignettes for more details about the entire process.
#'
#'
#' @section Algorithms:
#'   The model estimating functions in the \pkg{rstanarm} package take an 
#'   \code{algorithm} argument that can be one of the following:
#' \enumerate{
#'   \item \code{algorithm = "sampling"} Uses Markov Chain Monte Carlo
#'     (MCMC) --- in particular, Hamiltonian Monte Carlo (HMC) with a 
#'     tuned but diagonal mass matrix --- to draw from the posterior
#'     distribution of the parameters. See \code{\link[rstan]{sampling}}
#'     for more details. \strong{This is the default and the recommended
#'     algorithm for statistical inference.}
#'  \item \code{algorithm = "optimizing"} Finds the posterior mode using
#'    the LBGFS algorithm. See \code{\link[rstan]{optimizing}} for more details.
#'    If the priors are unspecified, then this is equivalent to maximum
#'    likelihood, in which case there is no great reason to use the functions in
#'    the \pkg{rstanarm} package over the emulated functions in other packages. 
#'    However, if priors are specified, then the estimates are penalized maximum
#'    likelihood estimates, which may have some redeeming value. 
#'    
#'    Optimization is known to perform poorly for mixed models and thus is
#'    not currently enabled for models fit with the \code{\link{stan_lmer}} 
#'    and \code{\link{stan_glmer}} functions (for these models 
#'    \code{algorithm="sampling"} should be used).
#' }
#' 
#' (Additional algorithms --- e.g., mean-field and full-rank variational
#' approximations to full Bayesian inference --- are in development and will be
#' available in future releases of the package.)
#' 
#' See \code{\link{priors}} for an overview of the various choices the
#' user can make for prior distributions. 
#' 
#' 
#' @section Modeling functions:
#' The model estimating functions are described in greater detail in their 
#' individual help pages and the vignettes. Here we provide a brief overview:
#' 
#' \describe{
#'  \item{\code{\link[=stan_lm]{stan_lm and stan_aov}}}{
#'   Similar to \code{\link[stats]{lm}} or \code{\link[stats]{aov}} but with 
#'   regularizing priors on the model parameters that are driven by prior 
#'   beliefs about \eqn{R^2}, the proportion of variance in the outcome 
#'   attributable to the predictors in a linear model.
#'  }
#'  \item{\code{\link{stan_glm}}}{
#'   Similar to \code{\link[stats]{glm}} but with Gaussian, Student t, Cauchy 
#'   or hierarhical prior distributions for the coefficients and a half-Cauchy
#'   prior for any nuisance parameter in a Generalized Linear Model (GLM) for
#'   various outcomes that are characterized by a \code{\link[stats]{family}}
#'   object. It is also possible to estimate a negative bionomial model in a 
#'   similar way to the \code{\link[MASS]{glm.nb}} function in the \pkg{MASS} 
#'   package.
#'  }
#'  \item{\code{\link[=stan_glmer]{stan_glmer and stan_lmer}}}{
#'   Similar to \code{\link[lme4]{glmer}} and \code{\link[lme4]{lmer}} in the
#'   \pkg{lme4} package, which augments GLMs to have group-specific terms that
#'   deviate from the common coefficients according to a mean-zero multivariate
#'   normal distribution with a highly-structured but unknown covariance matrix
#'   that itself has an innovative prior distribution. MCMC provides more 
#'   appropriate estimates of uncertainty of such models and provides posterior 
#'   distributions for each of the group-specific parameters.
#'  }
#'  \item{\code{\link{stan_polr}}}{
#'   Similar to \code{\link[MASS]{polr}} in the \pkg{MASS} package (in that it
#'   models an ordinal response in a similar way to \code{\link{stan_glm}}) but
#'   also specifies a prior on the unknown cutpoints. Can also be used to model
#'   binary outcomes.
#'  }
#' }
#'   
#' @seealso \code{\link{stanreg-objects}} and \code{\link{stanreg-methods}} for 
#'   details on the fitted model objects returned by the modeling functions.
#'   
#'   \code{\link{priors}} for an overview of the various choices the user can
#'   make for prior distributions.
#'   
#'   \code{\link{plots}} for the various plotting functions that can be used
#'   to explore fitted models.
#'   
#'   \url{http://mc-stan.org/} for more information on the Stan C++ package used
#'   by \pkg{rstanarm} to estimate the models.
#'   
NULL
