#' @param QR A logical scalar defaulting to \code{FALSE}, but if \code{TRUE}
#'   applies a scaled \code{\link{qr}} decomposition to the design matrix. The
#'   transformation does not change the likelihood of the data but is
#'   recommended for computational reasons when there are multiple predictors.
#'   See the \link{QR-argument} documentation page for details on how
#'   \pkg{rstanarm} does the transformation and important information about how
#'   to interpret the prior distributions of the model parameters when using
#'   \code{QR=TRUE}.
