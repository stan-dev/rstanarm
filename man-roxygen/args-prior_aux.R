#' @param prior_aux The prior distribution for the "auxiliary" parameter (if
#'   applicable). The "auxiliary" parameter refers to a different parameter 
#'   depending on the \code{family}. For Gaussian models \code{prior_aux} 
#'   controls \code{"sigma"}, the error 
#'   standard deviation. For negative binomial models \code{prior_aux} controls 
#'   \code{"reciprocal_dispersion"}, which is similar to the 
#'   \code{"size"} parameter of \code{\link[stats]{rnbinom}}:
#'   smaller values of \code{"reciprocal_dispersion"} correspond to 
#'   greater dispersion. For gamma models \code{prior_aux} sets the prior on 
#'   to the \code{"shape"} parameter (see e.g., 
#'   \code{\link[stats]{rgamma}}), and for inverse-Gaussian models it is the 
#'   so-called \code{"lambda"} parameter (which is essentially the reciprocal of
#'   a scale parameter). Binomial and Poisson models do not have auxiliary 
#'   parameters. 
#'   
#'   \code{prior_aux} can be a call to \code{exponential} to 
#'   use an exponential distribution, or \code{normal}, \code{student_t} or 
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy 
#'   prior. See \code{\link{priors}} for details on these functions. To omit a 
#'   prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_aux} to \code{NULL}.
