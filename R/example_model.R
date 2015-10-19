#' Example model
#' 
#' A pre-fit model for use in examples. 
#' 
#' @name example_model
#' @format A \code{\link[=stanreg-objects]{stanreg}} object containing the
#'   output from fitting the following model:
#'   
#' \code{stan_glmer(cbind(incidence, size - incidence) ~ size + period + (1 |
#' herd), data = lme4::cbpp, family = binomial)}
#' 
#' @seealso \code{\link[lme4]{cbpp}} for a description of the data.
#' 
NULL