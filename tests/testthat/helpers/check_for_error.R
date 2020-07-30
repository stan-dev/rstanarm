# These tests just make sure that posterior_predict doesn't throw errors and
# that result has correct dimensions
check_for_error <- function(fit, data = NULL, offset = NULL) {
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
