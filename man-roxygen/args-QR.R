#' @param QR A logical scalar (defaulting to \code{FALSE}) but if \code{TRUE}
#'   applies a scaled \code{\link{qr}} decomposition to the design matrix, 
#'   \eqn{X = Q* R*}{X = Q^\ast R^\ast}, where 
#'   \eqn{Q* = Q (n-1)^0.5}{Q^\ast = Q \sqrt{n-1}} and
#'   \eqn{R* = (n-1)^-0.5 R}{R^\ast = \frac{1}{\sqrt{n-1}} R}. The coefficients
#'   relative to \eqn{Q*}{Q^\ast} are obtained and then premultiplied by the
#'   inverse of \eqn{R*}{R^{\ast}} to obtain coefficients relative to the
#'   original predictors, \eqn{X}. These transformations do not change the 
#'   likelihood of the data but are recommended for computational reasons when 
#'   there are multiple predictors but you do not have an informative prior on
#'   their coefficients.
