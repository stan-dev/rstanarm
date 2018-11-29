# Recover parameter estimates and return a list with consistent
# parameter names for comparing stan_surv and coxph estimates
#
# @param mod The fitted survival model. Likely to be a model estimated 
#   using either coxph or stan_surv.
#
recover_pars <- function(mod) {

  cl <- class(mod)[1L]
  
  fixef_pars <- switch(cl,
                       coxph    = mod$coefficients,
                       survreg  = mod$coefficients,
                       stansurv = fixef(mod),
                       NULL)
  
  ret <- Filter(function(x) !is.null(x), list(fixef = fixef_pars))
  
  return(ret)
}
