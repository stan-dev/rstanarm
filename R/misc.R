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
  if (!is.null(x) & any(x <= 0)) {
    nm <- deparse(substitute(x))
    stop(paste(nm, "should be positive"))
  }
}

set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link))
  if (is.null(scale))
    scale <- default
  if (link == "probit")
    scale * dnorm(0) / dlogis(0)
  else 
    scale
}
