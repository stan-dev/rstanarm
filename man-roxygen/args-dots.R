#' @param ... Further arguments passed to the function in the \pkg{rstan} 
#'   package (\code{\link[rstan:stanmodel-method-sampling]{sampling}}, 
#'   \code{\link[rstan:stanmodel-method-vb]{vb}}, or 
#'   \code{\link[rstan:stanmodel-method-optimizing]{optimizing}}), 
#'   corresponding to the estimation method named by \code{algorithm}. For example, 
#'   if \code{algorithm} is \code{"sampling"} it is possible to specify \code{iter}, 
#'   \code{chains}, \code{cores}, and other MCMC controls.  
#'   
#'   Another useful argument that can be passed to \pkg{rstan} via \code{...} is
#'   \code{refresh}, which specifies how often to print updates when sampling
#'   (i.e., show the progress every \code{refresh} iterations). \code{refresh=0}
#'   turns off the iteration updates.
