#' @param QR A logical scalar (defaulting to \code{FALSE}) but if \code{TRUE}
#'   applies a scaled \code{\link{qr}} decomposition to the design matrix, 
#'   \eqn{X = Q^\ast R^\ast}{X = Q* R*}, where 
#'   \eqn{Q^\ast = Q \sqrt{n-1}}{Q* = Q (n-1)^0.5} and
#'   \eqn{R^\ast = \frac{1}{\sqrt{n-1}} R}{R* = (n-1)^(-0.5) R}. The coefficients
#'   relative to \eqn{Q^\ast}{Q*} are obtained and then premultiplied by the
#'   inverse of \eqn{R^{\ast}}{R*} to obtain coefficients relative to the
#'   original predictors, \eqn{X}. These transformations do not change the 
#'   likelihood of the data but are recommended for computational reasons when 
#'   there are multiple predictors. However, because the coefficients relative
#'   to \eqn{Q^\ast}{Q*} are not very interpretable it is hard to specify an 
#'   informative prior. Setting \code{QR=TRUE} is therefore only recommended 
#'   if you do not have an informative prior for the regression coefficients.
