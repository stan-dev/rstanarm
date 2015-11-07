#' Family function for Negative Binomial GLMs
#' 
#' Specifies the information required to fit a Negative Binomial GLM in a 
#' similar way to \code{\link[MASS]{negative.binomial}}. However, here the 
#' overdispersion parameter \code{theta} is not specified by the user and always
#' estimated. A call to this function can be passed to the \code{family}
#' argument of \code{\link{stan_glm}} or \code{\link{stan_glmer}} to estimate a
#' Negative Binomial model. Alternatively, the \code{\link{stan_glm.nb}} and 
#' \code{\link{stan_glmer.nb}} wrapper functions may be used, which call 
#' \code{neg_binomial_2} internally.
#' 
#' @export
#' @param link The same as for \code{\link{poisson}}, typically a character
#'   vector of length one among \code{"log"}, \code{"identity"}, and
#'   \code{"sqrt"}.
#' @return An object of class \code{\link[stats]{family}} very similar to
#'   that of \code{\link[stats]{poisson}} but with a different family name.
#' @examples
#' \dontrun{
#' # Example usage with stan_glm
#' options(mc.cores = 4)
#' SEED <- 12345
#' 
#' # can specify family = neg_binomial_2
#' fit <- stan_glm(Days ~ Sex/(Age + Eth*Lrn), data = MASS::quine, 
#'                 family = neg_binomial_2, seed = SEED) 
#'                 
#' # or, equivalently, use stan_glm.nb 
#' stan_glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = MASS::quine, 
#'             seed = SEED) 
#' }
#'
neg_binomial_2 <- function(link = "log") {
  out <- poisson(link)
  out$family <- "neg_binomial_2"
  out$variance <- function(mu, theta) mu + mu^2 / theta
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}
