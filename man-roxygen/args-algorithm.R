#' @param algorithm Character string (possibly abbreviated) indicating the 
#'   estimation approach to use. Can be \code{"sampling"} for MCMC (the
#'   default), \code{"optimizing"} for optimization, \code{"meanfield"} for
#'   variational inference with independent normal distributions, or
#'   \code{"fullrank"} for variational inference with a multivariate normal
#'   distribution. See \code{\link{rstanarm-package}} for more details on the
#'   estimation algorithms. NOTE: not all fitting functions support all four
#'   algorithms.
