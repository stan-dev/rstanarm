#' @param ... Further arguments passed to the function in the \pkg{rstan} 
#'   package (\code{\link[rstan]{sampling}} or \code{\link[rstan]{optimizing}},
#'   corresponding to the estimation method named by \code{algorithm}. For 
#'   example, if \code{algorithm='sampling'} we might specify \code{iter}, 
#'   \code{chains}, \code{cores}, etc.).
