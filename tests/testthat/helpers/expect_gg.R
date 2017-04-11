expect_gg <- function(x, info = NULL, label = NULL) {
  testthat::expect_is(x, "ggplot", info = info, label = label)
  invisible(ggplot2::ggplot_build(x))
}
