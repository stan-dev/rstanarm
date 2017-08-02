#' @param prior The prior distribution for the regression coefficients. 
#'   \code{prior} should be a call to one of the various functions provided by 
#'   \pkg{rstanarm} for specifying priors. The subset of these functions that 
#'   can be used for the prior on the coefficients can be grouped into several 
#'   "families":
#'   
#'   \tabular{ll}{
#'     \strong{Family} \tab \strong{Functions} \cr 
#'     \emph{Student t family} \tab \code{normal}, \code{student_t}, \code{cauchy} \cr 
#'     \emph{Hierarchical shrinkage family} \tab \code{hs}, \code{hs_plus} \cr 
#'     \emph{Laplace family} \tab \code{laplace}, \code{lasso} \cr
#'     \emph{Product normal family} \tab \code{product_normal} \cr
#'   }
#'   
#'   See the \link[=priors]{priors help page} for details on the families and 
#'   how to specify the arguments for all of the functions in the table above.
#'   To omit a prior ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior} can be set to \code{NULL}, although this is rarely a good
#'   idea.
#'   
#'   \strong{Note:} Unless \code{QR=TRUE}, if \code{prior} is from the Student t
#'   family or Laplace family, and if the \code{autoscale} argument to the 
#'   function used to specify the prior (e.g. \code{\link{normal}}) is left at 
#'   its default and recommended value of \code{TRUE}, then the default or 
#'   user-specified prior scale(s) may be adjusted internally based on the scales
#'   of the predictors. See the \link[=priors]{priors help page} for details on
#'   the rescaling and the \code{\link{prior_summary}} function for a summary of
#'   the priors used for a particular model.
#'   
