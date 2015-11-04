#' @param adapt_delta The target average proposal acceptance probability for 
#'   adaptation. The step size used by the numerical integrator is a function of
#'   \code{adapt_delta} in that increasing \code{adapt_delta} will result in a
#'   smaller step size. In general you should not need to change
#'   \code{adapt_delta} unless you see a warning message about divergent
#'   transitions, in which case you can increase \code{adapt_delta} from the
#'   default of 0.95 to a value closer to 1. Increasing \code{adapt_delta} will
#'   typically result in a slower sampler, but it will always lead to a more
#'   robust sampler.
