#' @param sparse A logical scalar (defaulting to \code{FALSE}) indicating
#'   whether to use a sparse representation of the design (X) matrix. 
#'   Setting this to \code{TRUE} will likely be twice as slow, even if the
#'   design matrix has a considerable number of zeros, but it may allow the
#'   model to be estimated when the computer has too little RAM to
#'   utilize a dense design matrix. If \code{TRUE}, the the design matrix
#'   is not centered (since that would destroy the sparsity) and it is
#'   not possible to specify both \code{QR = TRUE} and \code{sparse = TRUE}.
