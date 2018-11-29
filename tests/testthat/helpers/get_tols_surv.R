# Use the standard errors from a fitted 'comparison model' to obtain 
# the tolerance for each parameter in the joint model
# Obtain parameter specific tolerances that can be used to assess the 
# accuracy of parameter estimates in stan_jm models. The tolerances
# are calculated by taking the SE/SD for the parameter estimate in a 
# "gold standard" model and multiplying this by the relevant element 
# in the 'tolscales' argument.
#
# @param mod The "gold standard" longitudinal model. Likely to be
#   a model estimated using coxph.
# @param toscales A named list with elements 'hr_fixef' and 'tde_fixef'.
#
get_tols <- function(mod, tolscales) {
  
  cl <- class(mod)[1L]
  
  if (cl %in% c("coxph", "survreg")) {
    fixef_ses  <- sqrt(diag(mod$var))[1:length(mod$coefficients)]
    fixef_tols <- tolscales$hr_fixef * fixef_ses
    names(fixef_tols) <- names(mod$coefficients)
  }
   
  if ("(Intercept)" %in% names(fixef_tols))
    fixef_tols[["(Intercept)"]] <- 2 * fixef_tols[["(Intercept)"]]
  
  ret <- Filter(function(x) !is.null(x), list(fixef = fixef_tols))
  
  return(ret)
}
