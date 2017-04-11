expect_equivalent_loo <- function(fit) {
  l <- suppressWarnings(loo(fit))
  w <- suppressWarnings(waic(fit))
  expect_s3_class(l, "loo")
  expect_s3_class(w, "loo")
  expect_s3_class(w, "waic")
  
  att_names <- c("names", "log_lik_dim", "class", "name", "discrete", "yhash")
  expect_named(attributes(l), att_names)
  expect_named(attributes(w), att_names)
  
  discrete <- attr(l, "discrete")
  expect_true(!is.na(discrete) && is.logical(discrete))
  
  expect_equivalent(l, suppressWarnings(loo(log_lik(fit))))
  expect_equivalent(w, suppressWarnings(waic(log_lik(fit))))
}
