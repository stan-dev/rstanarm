#' Fitted model objects 
#' 
#' @name stanreg-objects 
#' 
#' @description The model-fitting functions in \pkg{rstanarm} return an object 
#'   of class 'stanreg', which is a list containing at a minimum the components 
#'   listed below. Each stanreg object will also have additional components and
#'   at least one additional class depending on the model. The additional
#'   class(es) can be 'lm', 'glm', 'polr', 'lmerMod', or 'aov'.
#'   
#' @section stanreg objects:   
#' \describe{
#'   \item{\code{coefficients}}{Point estimates. For MCMC these are posterior
#'   medians.}
#'   \item{\code{ses}}{Standard errors based on median absolute deviation
#'   (\code{\link[stats]{mad}}).}
#'   \item{\code{residuals}}{Residuals of type \code{'response'}. See
#'   \code{\link{residuals.stanreg}}.}
#'   \item{\code{fitted.values}}{Fitted mean values. For GLMs the linear
#'   predictors are transformed by the inverse link function.}
#'   \item{\code{linear.predictors}}{Linear fit on the link scale. For linear
#'   models this is the same as \code{fitted.values}.}
#'   \item{\code{covmat}}{Variance-covariance matrix for the coefficients. For
#'   MCMC this is estimated from the posterior draws.}
#'   \item{\code{model,x,y}}{If requested, the the model frame, model matrix and
#'   \code{y} vector used, respectively.}
#'   \item{\code{family}}{The \code{\link[stats]{family}} object used.}
#'   \item{\code{df.residual}}{the residual degrees of freedom}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{formula}}{The model \code{\link[stats]{formula}}.}
#'   \item{\code{terms}}{The \code{\link[stats]{terms}} object used.}
#'   \item{\code{data,offset,weights}}{The \code{data}, \code{offset}, and
#'   \code{weights} arguments.}
#'   \item{\code{algorithm}}{The estimation method used (e.g. "sampling").}
#'   \item{\code{prior.info}}{A list with information about the prior distributions
#'   used.}
#'   \item{\code{stanfit,stan_summary}}{The object of
#'   \code{\link[rstan]{stanfit-class}} returned by RStan and a matrix of
#'   various summary statistics from the stanfit object.}
#'   }
#'   
#' @seealso \code{\link{stanreg-methods}} 
#'   
NULL