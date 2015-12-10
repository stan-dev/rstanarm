#' Extract posterior sample
#' 
#' For models fit using MCMC (\code{algorithm='sampling'}), the posterior sample
#' is an \eqn{S} by \eqn{P} matrix, where \eqn{S} is the size of the sample --- 
#' the number of post-warmup draws from the posterior distribution of the model 
#' parameters --- and \eqn{P} is the number of parameters. If using optimization
#' (\code{algorithm='optimizing'}) or variational inference 
#' (\code{algorithm='meanfield'} or \code{algorithm='fullrank'}), there is no 
#' posterior sample but rather a \eqn{1000} by \eqn{P} matrix of draws from
#' either the asymptotic multivariate Gaussian sampling distribution of the
#' parameters or the variational approximation to the posterior distribution.
#' 
#' @method as.matrix stanreg
#' @export
#' 
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @param ... Optional arguments passed to \code{\link[rstan]{extract}}
#'   (\pkg{rstan}). Currently, only \code{pars} is supported, and only for
#'   models fit by MCMC.
#'   
#' @return A matrix or data frame, the dimensions of which depend on the model
#'   and estimation algorithm as described above.
#' 
#' @seealso \code{\link{stanreg-methods}}
#' 
#' @examples
#' # Extract posterior sample after MCMC
#' draws <- as.matrix(example_model)
#' 
#' # For example, we can see that the median of the draws for the intercept 
#' # is the same as the point estimate used for the intercept
#' print(median(draws[, "(Intercept)"]))
#' print(example_model$coefficients["(Intercept)"])
#' 
#' \dontrun{
#' # Extract draws from asymptotic Gaussian sampling distribution 
#' # after optimization
#' fit <- stan_glm(mpg ~ wt, data = mtcars, algorithm = "optimizing")
#' draws <- as.data.frame(fit)
#' print(colnames(draws))
#' 
#' # Extract draws from variational approximation to the posterior distribution
#' fit2 <- update(fit, algorithm = "meanfield")
#' draws <- as.data.frame(fit2)
#' print(colnames(draws))
#' }
#' 
as.matrix.stanreg <- function(x, ...) {
  msg <- "No draws found."
  if (used.optimizing(x)) {
    if ("pars" %in% names(list(...))) 
      message("'pars' argument ignored for models with algorithm='optimizing'.")
    out <- x$asymptotic_sampling_dist 
    if (is.null(out)) stop(msg)
    else {
      dispersion <- c("sigma", "scale", "lambda", "overdispersion")
      keep <- c(names(coef(x)), # return with coefficients first
                dispersion[which(dispersion %in% colnames(out))])
      return(out[, keep, drop = FALSE])
    }
  }
  sf <- x$stanfit
  if (sf@mode != 0) stop(msg)
  posterior <- rstan::extract(sf, permuted = FALSE, inc_warmup = FALSE, ...) 
  out <- apply(posterior, 3L, FUN = function(y) y)
  if (is(x, "lmerMod")) out <- unpad_reTrms(out, columns = TRUE)
  if (!"pars" %in% names(list(...))) {
    # remove mean_PPD and log-posterior unless user specified 'pars'
    sel <- !grepl("mean_PPD|log-posterior", colnames(out))
    out <- out[, sel, drop = FALSE]
  }
  return(out)
}

#' @rdname as.matrix.stanreg
#' @export
#' 
as.data.frame.stanreg <- function(x, ...) {
  as.data.frame(as.matrix.stanreg(x, ...))
}
