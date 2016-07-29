#' @param kappa_mean A positive scalar indicating the expectation of the concentration
#'   parameter in a von Mises-Fisher distribution where the random variable in question
#'   is each group's unit vector underlying its regression coefficients. See the vignette
#'   for \code{\link{stan_lm}} for more details. Defaults to \eqn{1} but is ignored if
#'   there is only one group. A higher value produces more shrinkage toward the common unit
#'   vector while a lower value makes the groups closer to independent.
