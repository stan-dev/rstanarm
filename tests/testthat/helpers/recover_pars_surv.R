# Recover parameter estimates and return a list with consistent
# parameter names for comparing stan_surv and coxph estimates
#
# @param mod The fitted survival model. Likely to be a model estimated 
#   using either coxph or stan_surv.
#
recover_pars <- function(mod) {

  cl <- class(mod)[1L]
  
  fixef_pars <- switch(cl,
                       "coxph"    = mod$coefficients,
                       "survreg"  = mod$coefficients,
                       "stansurv" = fixef(mod),
                       NULL)
  
  if (cl == "stansurv") {
    sel <- grep(":tve-[a-z][a-z]-coef[0-9]*$", names(fixef_pars))
    # replace stansurv tve names with coxph tt names
    if (length(sel)) {
      nms <- names(fixef_pars)[sel]
      nms <- gsub(":tve-[a-z][a-z]-coef[0-9]*$", "", nms)
      nms <- paste0("tt(", nms, ")")
      names(fixef_pars)[sel] <- nms 
    }
  }
  
  ret <- Filter(function(x) !is.null(x), list(fixef = fixef_pars))
  
  return(ret)
}
