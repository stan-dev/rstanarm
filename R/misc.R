na_replace <- function(x, replacement) {
  # if x is NA return replacement, else return x itself
  if (is.na(x)) 
    replacement 
  else 
    x
}

maybe_broadcast <- function(x, n) {
  # if x has length 1 replicate it n times, else return x itself
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

validate_parameter_value <- function(x) {
  if (!is.null(x) && x <= 0) {
    nm <- deparse(substitute(x))
    stop(paste(nm, "should be positive"))
  }
}
