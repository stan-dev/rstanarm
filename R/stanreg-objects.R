#' Fitted model objects 
#' 
#' @name stanreg-objects
#' @details The model-fitting functions in \pkg{rstanarm} return an object of 
#'   class \code{'stanreg'}, which is a list containing at a minimum the
#'   following components:
#' 
#' \describe{
#'   \item{\code{coefficients}}{Point estimates. For MCMC these are posterior means.}
#'   \item{\code{residuals}}{Residuals. See \code{\link{residuals.stanreg}}.}
#'   \item{\code{fitted.values}}{Fitted mean values. For GLMs the linear predictors are transformed by the invserse link function.}
#'   \item{\code{linear.predictors}}{Linear fit on the link scale. For linear models this is the same as \code{fitted.values}.}
#'   \item{\code{covmat}}{Variance-covariance matrix for the coefficients. For MCMC this is estimated from the posterior draws.}
#'   \item{\code{y}}{If requested, the \code{y} vector used.}
#'   \item{\code{x}}{If requested, the model matrix.}
#'   \item{\code{model}}{If requested, the model frame.}
#'   \item{\code{family}}{The \code{\link[stats]{family}} object used.}
#'   \item{\code{df.residual}}{the residual degrees of freedom}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{formula}}{The formula supplied.}
#'   \item{\code{data}}{The \code{data} argument.}
#'   \item{\code{algorithm}}{The estimation method used. Either "optimizing" or "sampling".}
#'   \item{\code{prior.info}}{A list with information about the prior distributions
#'   used.}
#'   \item{\code{stanfit}}{The object of \code{\link[rstan]{stanfit-class}} returned by RStan.}
#'   }
#' 
#' \subsection{Additional classes}{
#'  Each 'stanreg' object will also have at least one additional class. This
#'  will be one of 'lm', 'glm', 'polr', 'lmerMod', or 'aov', depending on the
#'  model.
#'  }
#'   
#' @seealso \code{\link{stanreg-methods}} 
#'   
NULL