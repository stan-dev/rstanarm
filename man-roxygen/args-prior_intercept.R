#' @param prior_intercept The prior distribution for the intercept (after
#'   centering all predictors, see note below).
#' 
#'   The default prior is described in the vignette 
#'   \href{https://mc-stan.org/rstanarm/articles/priors.html}{\emph{Prior
#'   Distributions for rstanarm Models}}.
#'   If not using the default, \code{prior_intercept} can be a call to
#'   \code{normal}, \code{student_t} or \code{cauchy}. See the
#'   \link[=priors]{priors help page} for details on these functions. To omit a
#'   prior on the intercept ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior_intercept} can be set to \code{NULL}.
#'   
#'   \strong{Note:} If using a dense representation of the design matrix
#'   ---i.e., if the \code{sparse} argument is left at its default value of
#'   \code{FALSE}--- then the prior distribution for the intercept is set so it
#'   applies to the value \emph{when all predictors are centered} (you don't
#'   need to manually center them). This is explained further in
#'   [Prior Distributions for rstanarm Models](https://mc-stan.org/rstanarm/articles/priors.html)
#'   If you prefer to specify a prior on the intercept without the predictors
#'   being auto-centered, then you have to omit the intercept from the
#'   \code{\link[stats]{formula}} and include a column of ones as a predictor,
#'   in which case some element of \code{prior} specifies the prior on it,
#'   rather than \code{prior_intercept}. Regardless of how
#'   \code{prior_intercept} is specified, the reported \emph{estimates} of the
#'   intercept always correspond to a parameterization without centered
#'   predictors (i.e., same as in \code{glm}).
#'   
#'   
