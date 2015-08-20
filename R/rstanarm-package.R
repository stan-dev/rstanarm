#' Applied Regression Modeling via RStan
#'
#' @docType package
#' @name rstanarm
#' 
#' @import rstan
#' @import stats
#' 
#' @export loo
#' @export waic
#' 
#' @description An extension to the \code{\link[rstan]{rstan}} package enabling
#'   the most commonly applied regression models to be specified using customary
#'   R modeling syntax. The user can choose to estimate a model by optimization
#'   or perform full Bayesian inference via Markov chain Monte Carlo sampling.
#' 
#' @section \code{\link[=stan_lm]{stan_lm and stan_aov}}: 
#' 
#'   Linear modeling with regularizing priors on the model parameters that are
#'   driven by prior beliefs about \eqn{R^2}, the proportion of variance in the
#'   outcome attributable to the predictors. See \code{\link{priors}} for an
#'   explanation of this critical point.
#'   
#' @section \code{\link{stan_glm}}:
#'   
#'   Generalized linear modeling with Gaussian, Student t, or Cauchy prior
#'   distributions for the coefficients.
#'   
#' @section \code{\link[=stan_glmer]{stan_glmer and stan_lmer}}:
#'   
#'   Generalized linear modeling with group-specific terms with Gaussian,
#'   Student t, or Cauchy prior distributions for the coefficients and flexible
#'   priors for the unknown covariance matrices.
#'   
#' @section \code{\link{stan_polr}}:
#'   
#'   Ordinal regression modeling with Gaussian, Student t, or Cauchy prior
#'   distributions for the coefficients.
#'   
#' 
#' @seealso \code{\link{stanreg-objects}}, \code{\link{stanreg-methods}}       
NULL
