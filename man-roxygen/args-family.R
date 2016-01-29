#' @param family There are two ways that the family argument can be specified. 
#'   The first way is exactly the same as how families are specified when 
#'   calling \code{\link[<%= pkg %>]{<%= pkgfun %>}}, e.g., \code{family = 
#'   gaussian(link = "log")}. The second method, a call to 
#'   \code{\link{rstanarm_family}}, allows the user to override the default 
#'   values for various family-specific parameters. Also, in addition to the
#'   traditional \code{\link[=family]{families}}, \pkg{rstanarm} supports
#'   families \code{\link{neg_binomial_2}} for negative binomial models and
#'   \code{\link{t_family}} for linear models with t-distributed errors.
