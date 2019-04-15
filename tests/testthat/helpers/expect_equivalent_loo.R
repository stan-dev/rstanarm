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
