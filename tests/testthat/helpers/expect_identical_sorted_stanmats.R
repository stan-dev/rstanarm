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
