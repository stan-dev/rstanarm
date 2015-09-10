# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)

context("ORifNULL")

test_that("%ORifNULL% works", {
  `%ORifNULL%` <- rstanarm:::`%ORifNULL%`
  a <- list(NULL, NA, NaN, 1, "a", FALSE, mat.or.vec(5,5))
  b <- 1
  ans <- c(1, a[-1])
  for (j in seq_along(a)) {
    expect_identical(a[[j]] %ORifNULL% b, ans[[j]])
  }
})

