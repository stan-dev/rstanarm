# Recover parameter estimates and return a list with consistent
# parameter names for comparing stan_jm, stan_mvmer, stan_{g}lmer,
# {g}lmer, and coxph estimates
#
# @param modLong The fitted longitudinal model. Likely to be
#   a model estimated using either {g}lmer or stan_{g}lmer.
# @param modEvent The fitted event model. Likely to be a model
#   estimated using coxph.
# @param idvar The name of the ID variable. Used to extract the estimates
#   for group-specific parameters that correspond to the individual/patient.
#
recover_pars <- function(modLong, modEvent = NULL, idvar = "id") {
  
  if (is.null(modEvent))
    modEvent <- modLong

  if (class(modLong)[1] %in% c("stanreg", "lmerMod", "glmerMod")) {
    fixef_pars <- fixef(modLong) 
    ranef_pars <- ranef(modLong)[[idvar]]
  } else if (class(modLong)[1] == "stanmvreg") {
    fixef_pars <- fixef(modLong)$Long1
    ranef_pars <- ranef(modLong)$Long1[[idvar]]
  }
  
  if (class(modEvent)[1] == "coxph") {
    event_pars <- modEvent$coefficients
  } else if (class(modEvent)[1] == "stanmvreg") {
    event_pars <- fixef(modEvent)$Event
  }
  
  list(fixef = fixef_pars, ranef = ranef_pars, event = event_pars)
}
