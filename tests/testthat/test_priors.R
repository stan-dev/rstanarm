library(rstanarm)
SEED <- 12345
set.seed(SEED)
ITER <- 10L
CHAINS <- 2L
REFRESH <- 0

SW <- suppressWarnings


context("priors")

is.rstanarm_prior <- rstanarm:::is.rstanarm_prior
el_nms <- c("dist", "df", "location", "scale")
test_that("normal works", {
  prior <- normal()
  expect_true(is.rstanarm_prior(prior))
  expect_identical(prior, normal(0, NULL))
  expect_identical(names(prior), el_nms)
  
  prior <- normal(2, 3)
  expect_identical(prior$df, NA)
  expect_identical(prior$location, 2)
  expect_identical(prior$scale, 3)
  
  prior <- normal(c(1,2), 3)
  expect_identical(prior$location, c(1,2))
  expect_identical(prior$scale, 3)
  
  expect_error(normal(0, -1), "positive")
})

test_that("student_t works", {
  prior <- student_t()
  expect_true(is.rstanarm_prior(prior))
  expect_identical(prior, student_t(1, 0, NULL))
  expect_identical(names(prior), el_nms)
  
  prior <- student_t(3, 2, 1)
  expect_identical(prior$df, 3)
  expect_identical(prior$location, 2)
  expect_identical(prior$scale, 1)
  
  prior <- student_t(location = c(1,2), scale = 3)
  expect_identical(prior$df, 1)
  expect_identical(prior$location, c(1,2))
  expect_identical(prior$scale, 3)
  
  expect_error(student_t(0, -1), "positive")
  expect_error(student_t(0, 0, 1), "positive")
})

test_that("cauchy works", {
  prior <- cauchy()
  expect_true(is.rstanarm_prior(prior))
  expect_identical(prior, student_t(1, 0, NULL))
  expect_identical(names(prior), el_nms)
  
  prior <- cauchy(2, 1)
  expect_identical(prior$df, 1)
  expect_identical(prior$location, 2)
  expect_identical(prior$scale, 1)
  
  prior <- cauchy(location = c(1,2), scale = 3)
  expect_identical(prior$df, 1)
  expect_identical(prior$location, c(1,2))
  expect_identical(prior$scale, 3)
  
  expect_error(cauchy(0, -1), "positive")
})

test_that("hs works", {
  prior <- hs()
  expect_true(is.rstanarm_prior(prior))
  expect_identical(names(prior), el_nms)
  expect_identical(prior$dist, "hs")
  expect_identical(prior$df, 3)
  expect_identical(prior$location, 0)
  expect_identical(prior$scale, 1)
  
  prior <- hs(5)
  expect_identical(prior$df, 5)
  expect_identical(prior$location, 0)
  expect_identical(prior$scale, 1)

  expect_error(hs(0), "positive")
  expect_error(hs(-1), "positive")
})

test_that("hs_plus works", {
  prior <- hs_plus()
  expect_true(is.rstanarm_prior(prior))
  expect_identical(names(prior), el_nms)
  expect_identical(prior$dist, "hs_plus")
  expect_identical(prior$df, 3)
  expect_identical(prior$location, 0)
  expect_identical(prior$scale, 3) # scale used for second df parameter
  
  prior <- hs_plus(5, 2)
  expect_identical(prior$df, 5)
  expect_identical(prior$location, 0)
  expect_identical(prior$scale, 2)
  
  prior <- hs_plus(df2 = 2)
  expect_identical(prior$df, 3)
  expect_identical(prior$scale, 2)
  
  expect_error(hs_plus(0), "positive")
  expect_error(hs_plus(3, -1), "positive")
})

test_that("dirichlet works", {
  prior <- dirichlet()
  expect_true(is.rstanarm_prior(prior))
  expect_identical(prior$dist, "dirichlet")
  expect_identical(prior$concentration, 1)
  
  prior <- dirichlet(concentration = 5)
  expect_identical(prior$concentration, 5)
  
  prior <- dirichlet(concentration = 1:5)
  expect_identical(prior$concentration, 1:5)
  
  expect_error(dirichlet(0))
  expect_error(dirichlet(-1))
  expect_error(dirichlet(c(-1, 2)))
})


test_that("decov works", {
  prior <- decov()
  el_nms <- c("dist", "regularization", "concentration", "shape", "scale")
  expect_true(is.rstanarm_prior(prior))
  expect_identical(prior, decov(1,1,1,1))
  expect_identical(names(prior), el_nms)

  expect_error(decov(0))
  expect_error(decov(scale = 0))
  expect_error(decov(-1,1,1,1))
  expect_error(decov(2,0,2,1))
})

test_that("R2 works", {
  prior <- R2(location = 0.5)
  el_nms <- c("dist", "location", "what", "df", "scale")
  expect_true(is.rstanarm_prior(prior))
  expect_identical(names(prior), el_nms)
  expect_identical(prior, R2(location = 0.5, what = "mode"))
  
  prior <- R2(0.25, what = "mean")
  expect_identical(prior$location, 0.25)
  expect_identical(prior$what, "mean")
  
  prior <- R2(-1, what = "log")
  expect_identical(prior$location, -1)
  expect_identical(prior$what, "log")
  
  expect_error(R2(), "'location' must be numeric")
  expect_error(R2(0.5, what = "banana"), "should be one of")
  for (wh in c("log", "mode", "mean", "median")) 
    expect_error(R2(2, what = wh), "'location' must")
  for (wh in c("mode", "mean", "median")) 
    expect_error(R2(-1, what = wh), "'location' must")
})


context("priors (helpers)")

test_that("set_prior_scale works", {
  set_prior_scale <- rstanarm:::set_prior_scale
  expect_error(set_prior_scale("a", "b", "c"))
  expect_error(set_prior_scale(1, 1, 1))
  expect_equal(set_prior_scale(NULL, 1, "a"), 1)
  expect_equal(set_prior_scale(NULL, 1, "probit"), dnorm(0) / dlogis(0))
  expect_equal(set_prior_scale(2, 1, "a"), 2)
  expect_equal(set_prior_scale(2, 1, "probit"), 2 * dnorm(0) / dlogis(0))
})
test_that("validate_parameter_value works", {
  validate_parameter_value <- rstanarm:::validate_parameter_value
  expect_error(validate_parameter_value(-1), "should be positive")
  expect_error(validate_parameter_value(0), "should be positive")
  expect_error(validate_parameter_value("a"), "should be NULL or numeric")
  expect_error(validate_parameter_value(NA), "should be NULL or numeric")
  expect_true(validate_parameter_value(NULL))
  expect_true(validate_parameter_value(.01))
  expect_true(validate_parameter_value(.Machine$double.xmax))
})

test_that("validate_glm_prior works", {
  validate_glm_prior <- rstanarm:::validate_glm_prior
  
  # prior on coefficients
  val <- validate_glm_prior(normal(0,1), prior_for = "coef", link = "identity", ncoef = 5)
  ans <- list(dist = 1, scale = rep(1, 5), mean = array(rep(0, 5)), df = array(rep(1, 5)))
  expect_equal(val, ans)
  
  val <- validate_glm_prior(student_t(3), prior_for = "coef", link = "probit", ncoef = 5)
  ans <- list(dist = 2, scale = NULL, mean = array(rep(0, 5)), df = array(rep(3, 5)))
  ans$scale <- rep(rstanarm:::set_prior_scale(ans$scale, default = 2.5, link = "probit"), 5)
  expect_equal(val, ans)
  
  val <- validate_glm_prior(hs_plus(4,4), prior_for = "coef", link = "log", ncoef = 10)
  ans <- list(dist = 4, scale = rep(4, 10), mean = array(rep(0, 10)), df = array(rep(4, 10)))
  expect_equal(val, ans)
  
  val <- validate_glm_prior(normal(1:4, 2), prior_for = "coef", link = "log", ncoef = 4)
  ans <- list(dist = 1, scale = rep(2, 4), mean = array(1:4), df = array(rep(1, 4)))
  expect_equal(val, ans)
  
  # prior on intercept
  val <- validate_glm_prior(normal(), prior_for = "intercept", link = "log")
  ans <- list(dist = 1, scale = 10, mean = 0, df = 1)
  expect_equal(val, ans)
  
  val <- validate_glm_prior(normal(-1, 3), prior_for = "intercept", link = "logit")
  ans <- list(dist = 1, scale = 3, mean = -1, df = 1)
  expect_equal(val, ans)
  
  val <- validate_glm_prior(cauchy(0, 2.5), prior_for = "intercept", link = "probit")
  ans <- list(dist = 2, scale = 2.5, mean = 0, df = 1)
  ans$scale <- rstanarm:::set_prior_scale(ans$scale, default = 10, link = "probit")
  expect_equal(val, ans)
  
  val <- validate_glm_prior(student_t(7), prior_for = "intercept", link = "probit")
  ans <- list(dist = 2, scale = NULL, mean = 0, df = 7)
  ans$scale <- rstanarm:::set_prior_scale(ans$scale, default = 10, link = "probit")
  expect_equal(val, ans)
  
  expect_error(validate_glm_prior(hs(3), prior_for = "intercept", link = "identity"), 
               regexp = "distribution for the intercept should be one of")
  expect_error(validate_glm_prior(R2(0.5), prior_for = "coef", link = "identity", ncoef = 2), 
               regexp = "distribution for the coefficients should be one of")
  expect_error(validate_glm_prior(letters, prior_for = "intercept", link = "identity"), 
               regexp = "should be a named list")
  expect_error(validate_glm_prior(letters, prior_for = "coef", link = "identity", ncoef = 2), 
               regexp = "should be a named list")
})

