expect_ppd <- function(x) {
  expect_true(inherits(x, "ppd") || is.matrix(x))
}
