#' @param sparse A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to use a sparse representation of the design (X) matrix. 
#'   If \code{TRUE}, the the design matrix is not centered (since that would 
#'   destroy the sparsity) and likewise it is not possible to specify both 
#'   \code{QR = TRUE} and \code{sparse = TRUE}. Depending on how many zeros
#'   there are in the design matrix, setting \code{sparse = TRUE} may make
#'   the code run faster and can consume much less RAM. 
