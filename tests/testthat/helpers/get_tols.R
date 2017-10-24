# Use the standard errors from a fitted 'comparison model' to obtain 
# the tolerance for each parameter in the joint model
# Obtain parameter specific tolerances that can be used to assess the 
# accuracy of parameter estimates in stan_jm models. The tolerances
# are calculated by taking the SE/SD for the parameter estimate in a 
# "gold standard" model and multiplying this by the relevant element 
# in the 'tolscales' argument.
#
# @param modLong The "gold standard" longitudinal model. Likely to be
#   a model estimated using either {g}lmer or stan_{g}lmer.
# @param modEvent The "gold standard" event model. Likely to be a model
#   estimated using coxph.
# @param toscales A named list with elements $lmer_fixef, $lmer_ranef,
#   $glmer_fixef, $glmer_ranef, $event.
# @param idvar The name of the ID variable. Used to extract the SDs for
#   group-specific terms that correspond to the individual/patient.
#
get_tols <- function(modLong, modEvent = NULL, tolscales, idvar = "id") {
  
  if (is.null(modEvent))
    modEvent <- modLong # if modLong is already a joint model
  
  if (class(modLong)[1] == "stanreg") {
    fixef_nms <- names(fixef(modLong))
    fixef_ses <- modLong$ses[fixef_nms] 
    ranef_sds <- attr(VarCorr(modLong)[[idvar]], "stddev")
    if (modLong$stan_function == "stan_lmer") {
      fixef_tols <- tolscales$lmer_fixef * fixef_ses
      ranef_tols <- tolscales$lmer_ranef * ranef_sds
    } else if (modLong$stan_function == "stan_glmer") {
      if (modLong$family$family == "gaussian") {
        fixef_tols <- tolscales$lmer_fixef * fixef_ses
        ranef_tols <- tolscales$lmer_ranef * ranef_sds
      } else {
        fixef_tols <- tolscales$glmer_fixef * fixef_ses
        ranef_tols <- tolscales$glmer_ranef * ranef_sds
      }
    }
  } else if (class(modLong)[1] %in% c("lmerMod", "glmerMod")) {
    fixef_ses <- sqrt(diag(vcov(modLong)))
    ranef_sds <- attr(VarCorr(modLong)[[idvar]], "stddev")
    if (class(modLong)[1] == "lmerMod") {
      fixef_tols <- tolscales$lmer_fixef * fixef_ses
      ranef_tols <- tolscales$lmer_ranef * ranef_sds      
    } else if (class(modLong)[1] == "glmerMod") {
      fixef_tols <- tolscales$glmer_fixef * fixef_ses
      ranef_tols <- tolscales$glmer_ranef * ranef_sds
    }
  }
  if ("(Intercept)" %in% names(fixef_tols))
    fixef_tols[["(Intercept)"]] <- 2 * fixef_tols[["(Intercept)"]]
  if ("(Intercept)" %in% names(ranef_tols))
    ranef_tols[["(Intercept)"]] <- 2 * ranef_tols[["(Intercept)"]]
  
  if (class(modEvent)[1] == "coxph") {
    event_ses <- summary(modEvent)$coefficients[, "se(coef)"]
  } else event_ses <- NULL
  event_tols <- if (!is.null(event_ses))
    tolscales$event * event_ses else NULL
  if ("(Intercept)" %in% names(event_tols))
    event_tols[["(Intercept)"]] <- 2 * event_tols[["(Intercept)"]]
  
  ret <- Filter(
    function(x) !is.null(x),
    list(fixef = fixef_tols, ranef = ranef_tols, event = event_tols))
  return(ret)
}
