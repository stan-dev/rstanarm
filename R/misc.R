na_replace <- function(x, replacement) {
  # return x itself if not na, else return replacement
  if (is.na(x)) 
    replacement 
  else 
    x
}

maybe_broadcast <- function(x, n) {
  # if x has length 1 it is replicated n times, else x itself is returned
  if (length(x) == 1L) 
    rep(x, times = n)
  else 
    x
}

is.binfac <- function(x) {
  # test if x is a factor with 2 levels (binary factor)
  is.factor(x) && nlevels(x) == 2L
}

fac2bin <- function(x) {
  # convert factor with 2 levels to 0/1
  if (!is.binfac(x)) 
    stop("x should be a factor with 2 levels")
  z <- as.numeric(x)
  ifelse(z == max(z), 1L, 0L)
}

nlist <- function(...) {
  # named lists
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names)
    FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }
  out
}

validate_loc_scale_df <- function(location, scale, df) {
  testthat::expect_is(location, "numeric")
  if (!is.null(scale)) 
    testthat::expect_more_than(scale, 0)
  if (!missing(df)) 
    testthat::expect_more_than(df, 0)
}
