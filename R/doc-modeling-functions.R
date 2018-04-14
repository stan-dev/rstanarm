#' Modeling functions available in \pkg{rstanarm}
#' 
#' @name available-models
#' 
#' @section Modeling functions:
#' The model estimating functions are described in greater detail in their
#' individual help pages and vignettes. Here we provide a very brief
#' overview:
#'
#' \describe{
#'  \item{\code{\link{stan_lm}}, \code{stan_aov}, \code{stan_biglm}}{
#'   Similar to \code{\link[stats]{lm}} or \code{\link[stats]{aov}} but with
#'   novel regularizing priors on the model parameters that are driven by prior
#'   beliefs about \eqn{R^2}, the proportion of variance in the outcome
#'   attributable to the predictors in a linear model.
#'  }
#'  \item{\code{\link{stan_glm}}, \code{stan_glm.nb}}{
#'  Similar to \code{\link[stats]{glm}} but with various possible prior
#'  distributions for the coefficients and, if applicable, a prior distribution
#'  for any auxiliary parameter in a Generalized Linear Model (GLM) that is
#'  characterized by a \code{\link[stats]{family}} object (e.g. the shape
#'  parameter in Gamma models). It is also possible to estimate a negative
#'  binomial model in a similar way to the \code{\link[MASS]{glm.nb}} function
#'  in the \pkg{MASS} package.
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
#'  \item{\code{\link{stan_nlmer}}}{
#'   Similar to \code{\link[lme4]{nlmer}} in the \pkg{lme4} package for 
#'   nonlinear "mixed-effects" models, but the group-specific coefficients 
#'   have flexible priors on their unknown covariance matrices.
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
#'   models an ordinal response, but the Bayesian model also implies a prior
#'   distribution on the unknown cutpoints. Can also be used to model binary
#'   outcomes, possibly while estimating an unknown exponent governing the
#'   probability of success.
#'  }
#'  \item{\code{\link{stan_betareg}}}{
#'   Similar to \code{\link[betareg]{betareg}} in that it models an outcome that
#'   is a rate (proportion) but, rather than performing maximum likelihood
#'   estimation, full Bayesian estimation is performed by default, with
#'   customizable prior distributions for all parameters.
#'  }
#'  \item{\code{\link{stan_clogit}}}{
#'    Similar to \code{\link[survival]{clogit}} in that it models an binary outcome
#'    where the number of successes and failures is fixed within each stratum by
#'    the research design. There are some minor syntactical differences relative
#'    to \code{\link[survival]{clogit}} that allow \code{stan_clogit} to accept
#'    group-specific terms as in \code{\link{stan_glmer}}.
#'  }
#'  \item{\code{\link{stan_mvmer}}}{
#'    A multivariate form of \code{\link{stan_glmer}}, whereby the user can
#'    specify one or more submodels each consisting of a GLM with group-specific 
#'    terms. If more than one submodel is specified (i.e. there is more than one
#'    outcome variable) then a dependence is induced by assuming that the
#'    group-specific terms for each grouping factor are correlated across submodels. 
#'  }
#'  \item{\code{\link{stan_jm}}}{
#'    Estimates shared parameter joint models for longitudinal and time-to-event 
#'    (i.e. survival) data. The joint model can be univariate (i.e. one longitudinal 
#'    outcome) or multivariate (i.e. more than one longitudinal outcome). A variety 
#'    of parameterisations are available for linking the longitudinal and event 
#'    processes (i.e. a variety of association structures).      
#'  }
#' }
#' 
#' @seealso \url{http://mc-stan.org/rstanarm/}
#' 
NULL
