#' Extract posterior sample
#' 
#' For models fit using MCMC (\code{algorithm='sampling'}), the posterior sample
#' is an \eqn{S} by \eqn{P} matrix, where \eqn{S} is the size of the sample --- 
#' the number of post-warmup draws from the posterior distribution of the model 
#' parameters --- and \eqn{P} is the number of parameters. For models fit via
#' optimization (\code{algorithm='optimizing'}), there is no posterior sample
#' but rather a \eqn{1000} by \eqn{P} matrix of draws from the asymptotic 
#' multivariate Gaussian sampling distribution of the parameters.
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
#' @return A matrix, the dimensions of which depend on the model and estimation
#'   algorithm as described above.
#' 
#' @seealso \code{\link{stanreg-methods}}
#' 
#' @examples
#' draws <- as.matrix(example_model)
#' 
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
  stopifnot(used.sampling(x))
  sf <- x$stanfit
  if (sf@mode != 0) stop(msg)
  posterior <- rstan::extract(sf, permuted = FALSE, inc_warmup = FALSE, ...) 
  out <- apply(posterior, 3L, FUN = function(y) y)
  if (is(x, "lmerMod")) out <- unpad_reTrms(out, columns = TRUE)
  # remove mean_PPD and log-posterior unless user specified 'pars'
  if ("pars" %in% names(list(...))) return(out)
  else out[, !grepl("mean_PPD|log-posterior", colnames(out)), drop = FALSE]
}
