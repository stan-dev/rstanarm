#' @param adapt_delta Only relevant if \code{algorithm='sampling'}.
#'   \code{adapt_delta} is the target average proposal acceptance probability
#'   for adaptation. The step size used by the numerical integrator is a
#'   function of \code{adapt_delta} in that increasing \code{adapt_delta} will
#'   result in a smaller step size. The default value of \code{adapt_delta} is
#'   0.95, except when the prior for the regression coefficients is
#'   \code{\link{R2}}, \code{\link{hs}}, \code{\link{hs_plus}}, or
#'   \code{\link{student_t}} with \code{df <= 2}, in which case the default is
#'   0.99. In general you should not need to change \code{adapt_delta} unless
#'   you see a warning message about divergent transitions, in which case you
#'   can increase \code{adapt_delta} from the default to a value \emph{closer}
#'   to 1 (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc). Increasing 
#'   \code{adapt_delta} will typically result in a slower sampler, but it will 
#'   always lead to a more robust sampler.
