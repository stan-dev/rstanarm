#' @param ... Further arguments passed to the function in the \pkg{rstan} 
#'   package (\code{\link[rstan]{sampling}}, \code{\link[rstan]{vb}}, or 
#'   \code{\link[rstan]{optimizing}}), corresponding to the estimation method 
#'   named by \code{algorithm}. For example, if \code{algorithm} is
#'   \code{"sampling"} it is possibly to specify \code{iter}, \code{chains},
#'   \code{cores}, \code{refresh}, etc.
