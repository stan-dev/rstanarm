SW <- function(expr) utils::capture.output(suppressWarnings(expr))

run_example_model <- function() {
  o <- SW(
    fit <- stan_glmer(cbind(incidence, size - incidence) ~ size + period + (1|herd),
             data = lme4::cbpp, family = binomial, QR = TRUE,
             # this next line is only to keep the example small in size!
             chains = 2, cores = 1, seed = 12345, iter = 1000, refresh = 0)
  )
  fit
}

# These tests just make sure that posterior_predict doesn't throw errors and
# that result has correct dimensions
check_for_pp_errors <- function(fit, data = NULL, offset = NULL) {
  nsims <- nrow(as.data.frame(fit))
  mf <- if (!is.null(data)) 
    data else model.frame(fit)
  if (identical(deparse(substitute(fit)), "example_model"))
    mf <- lme4::cbpp
  
  expect_silent(yrep1 <- posterior_predict(fit))
  expect_silent(lin1 <- posterior_linpred(fit))
  expect_silent(suppressMessages(posterior_linpred(fit, transform = TRUE)))
  expect_equal(dim(yrep1), c(nsims, nobs(fit)))
  expect_equal(dim(lin1), c(nsims, nobs(fit)))
  
  expect_silent(yrep2 <- posterior_predict(fit, draws = 1))
  expect_equal(dim(yrep2), c(1, nobs(fit)))
  
  offs <- if (!is.null(offset)) offset[1] else offset
  expect_silent(yrep3 <- posterior_predict(fit, newdata = mf[1,], offset = offs))
  expect_silent(lin3 <- posterior_linpred(fit, newdata = mf[1,], offset = offs))
  expect_equal(dim(yrep3), c(nsims, 1))
  expect_equal(dim(lin3), c(nsims, 1))
  
  expect_silent(yrep4 <- posterior_predict(fit, draws = 2, newdata = mf[1,], offset = offs))
  expect_equal(dim(yrep4), c(2, 1))
  
  offs <- if (!is.null(offset)) offset[1:5] else offset
  expect_silent(yrep5 <- posterior_predict(fit, newdata = mf[1:5,], offset = offs))
  expect_silent(lin5 <- posterior_linpred(fit, newdata = mf[1:5,], offset = offs))
  expect_equal(dim(yrep5), c(nsims, 5))
  expect_equal(dim(lin5), c(nsims, 5))
  
  expect_silent(yrep6 <- posterior_predict(fit, draws = 3, newdata = mf[1:5,], offset = offs))
  expect_equal(dim(yrep6), c(3, 5))
  
  expect_error(posterior_predict(fit, draws = nsims + 1), 
               regexep = "posterior sample size is only")
}


expect_equivalent_loo <- function(fit) {
  LOO.CORES <- ifelse(.Platform$OS.type == "windows", 1, 2)
  l <- suppressWarnings(loo(fit, cores = LOO.CORES))
  w <- suppressWarnings(waic(fit))
  expect_s3_class(l, "psis_loo")
  expect_s3_class(l, "loo")
  expect_s3_class(w, "loo")
  expect_s3_class(w, "waic")
  
  att_names <- c("names", "dims", "class", "model_name", "discrete", "yhash", "formula")
  expect_named(attributes(l), att_names)
  expect_named(attributes(w), att_names)
  
  discrete <- attr(l, "discrete")
  expect_true(!is.na(discrete) && is.logical(discrete))
  
  if (fit$stan_function != "stan_clogit") {
    ll <- log_lik(fit)
    r_eff <- loo::relative_eff(exp(ll), chain_id = rstanarm:::chain_id_for_loo(fit))
    l2 <- suppressWarnings(loo(ll, r_eff = r_eff, cores = LOO.CORES))
    expect_equal(l$estimates, l2$estimates)
    expect_equivalent(w, suppressWarnings(waic(ll)))
  }
}

expect_gg <- function(x, info = NULL, label = NULL) {
  testthat::expect_is(x, "ggplot", info = info, label = label)
  invisible(ggplot2::ggplot_build(x))
}


# Make sure that the fitted Stan models x and y have identical MCMC samples
# after sorting the stanmat columns (ie. parameters) by name
expect_identical_sorted_stanmats <- function(x, y) {
  x_mat <- as.matrix(x)
  y_mat <- as.matrix(y)
  x_nms <- colnames(x_mat)
  y_nms <- colnames(y_mat)
  x_mat_sorted <- x_mat[, order(x_nms), drop = FALSE]
  y_mat_sorted <- y_mat[, order(y_nms), drop = FALSE]
  expect_identical(x_mat_sorted, y_mat_sorted)
}

expect_linpred_equal <- function(object, tol = 0.1) {
  linpred <- posterior_linpred(object)
  expect_equal(apply(linpred, 2, median), object$linear.predictors, 
               tolerance = tol, 
               check.attributes = FALSE)
}

expect_matrix <- function(x) expect_true(is.matrix(x))

expect_ppd <- function(x) {
  expect_true(inherits(x, "ppd") || is.matrix(x))
}

expect_stanreg <- function(x) expect_s3_class(x, "stanreg")
expect_stanmvreg  <- function(x) expect_s3_class(x, "stanmvreg")
expect_survfit <- function(x) expect_s3_class(x, "survfit.stanjm")

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
  } else if (class(modLong)[1] %in% c("stanjm", "stanmvreg")) {
    fixef_pars <- fixef(modLong)[[1L]]
    ranef_pars <- ranef(modLong)[[1L]][[idvar]]
  }
  
  if (class(modEvent)[1] == "coxph") {
    event_pars <- modEvent$coefficients
  } else if (class(modEvent)[1] %in% c("stanjm", "stanmvreg")) {
    event_pars <- fixef(modEvent)$Event
  } else event_pars <- NULL
  
  ret <- Filter(
    function(x) !is.null(x),
    list(fixef = fixef_pars, ranef = ranef_pars, event = event_pars))
  return(ret)
}

