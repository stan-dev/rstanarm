#' @param prior_smooth The prior distribution for the hyperparameters in GAMs,
#'   with lower values yielding less flexible smooth functions.
#'    
#'   \code{prior_smooth} can be a call to \code{exponential} to 
#'   use an exponential distribution, or \code{normal}, \code{student_t} or 
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy 
#'   prior. See \code{\link{priors}} for details on these functions. To omit a 
#'   prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_smooth} to \code{NULL}. The number of hyperparameters depends
#'   on the model specification but a scalar prior will be recylced as necessary
#'   to the appropriate length.
