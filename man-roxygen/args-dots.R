#' @param ... Further arguments passed to the function in the \pkg{rstan} 
#'   package (\code{\link[rstan:stanmodel-method-sampling]{sampling}}, 
#'   \code{\link[rstan:stanmodel-method-vb]{vb}}, or 
#'   \code{\link[rstan:stanmodel-method-optimizing]{optimizing}}), 
#'   corresponding to the estimation method named by \code{algorithm}. For example, 
#'   if \code{algorithm} is \code{"sampling"} it is possibly to specify \code{iter}, 
#'   \code{chains}, \code{cores}, \code{refresh}, etc.
