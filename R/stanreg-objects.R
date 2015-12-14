#' Fitted model objects 
#' 
#' @name stanreg-objects 
#' 
#' @description The model-fitting functions in \pkg{rstanarm} return an object 
#'   of class 'stanreg', which is a list containing at a minimum the components 
#'   listed below. Each stanreg object will also have additional classes (e.g. 
#'   'aov', 'glm', 'polr', etc.) and additional components depending on the 
#'   model and estimation algorithm.
#'   
#' @section stanreg objects:   
#' \describe{
#'   \item{\code{coefficients}}{
#'   Point estimates, as described in \code{\link{print.stanreg}}.
#'   }
#'   \item{\code{ses}}{
#'   Standard errors based on \code{\link[stats]{mad}}, as described in
#'   \code{\link{print.stanreg}}.
#'   }
#'   \item{\code{residuals}}{
#'   Residuals of type \code{'response'}.
#'   }
#'   \item{\code{fitted.values}}{
#'   Fitted mean values. For GLMs the linear predictors are transformed by the
#'   inverse link function.
#'   }
#'   \item{\code{linear.predictors}}{
#'   Linear fit on the link scale. For linear models this is the same as
#'   \code{fitted.values}.
#'   }
#'   \item{\code{covmat}}{
#'   Variance-covariance matrix for the coefficients based on draws from the
#'   posterior distribution, the variational approximation, or the asymptotic 
#'   sampling distribution, depending on the estimation algorithm.
#'   }
#'   \item{\code{model,x,y}}{
#'   If requested, the the model frame, model matrix and response variable used, 
#'   respectively.
#'   }
#'   \item{\code{family}}{
#'   The \code{\link[stats]{family}} object used.
#'   }
#'   \item{\code{df.residual}}{
#'   The residual degrees of freedom.
#'   }
#'   \item{\code{call}}{
#'   The matched call.
#'   }
#'   \item{\code{formula}}{
#'   The model \code{\link[stats]{formula}}.
#'   }
#'   \item{\code{data,offset,weights}}{
#'   The \code{data}, \code{offset}, and \code{weights} arguments.
#'   }
#'   \item{\code{algorithm}}{
#'   The estimation method used.
#'   }
#'   \item{\code{prior.info}}{
#'   A list with information about the prior distributions used.
#'   }
#'   \item{\code{stanfit,stan_summary}}{
#'   The object of \code{\link[rstan]{stanfit-class}} returned by RStan and a
#'   matrix of various summary statistics from the stanfit object.
#'   }
#'}
#'
#' @seealso \code{\link{stanreg-methods}} 
#'   
NULL