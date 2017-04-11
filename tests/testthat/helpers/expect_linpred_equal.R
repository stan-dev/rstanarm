expect_linpred_equal <- function(object, tol = 0.1) {
  linpred <- posterior_linpred(object)
  expect_equal(apply(linpred, 2, median), object$linear.predictors, 
               tolerance = tol, 
               check.attributes = FALSE)
}
