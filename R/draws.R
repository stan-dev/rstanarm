#' Create a \code{draws} object from a \code{stanreg} object
#'
#' Convert a \code{stanreg} object to a format supported by the
#' \pkg{\link[posterior:posterior-package]{posterior}} package. 
#'
#' @name stanreg-draws-formats
#' @aliases as_draws as_draws_matrix as_draws_array as_draws_df as_draws_rvars as_draws_list
#' 
#' @param x A \code{stanreg} object returned by one of the \pkg{rstanarm}
#'   modeling functions.
#' @param ... Arguments (e.g., \code{pars}, \code{regex_pars}) passed internally to
#'   \code{\link{as.matrix.stanreg}} or \code{as.array.stanreg}.
#'   
#' @details To subset iterations, chains, or draws, use
#'   \code{\link[posterior:subset_draws]{subset_draws}} after making the
#'   \code{draws} object. To subset variables use \code{...} to pass the \code{pars}
#'   and/or \code{regex_pars} arguments to \code{as.matrix.stanreg} or
#'   \code{as.array.stanreg} (these are called internally by
#'   \code{as_draws.stanreg}), or use
#'   \code{\link[posterior:subset_draws]{subset_draws}} after making the
#'   \code{draws} object.
#'   
#' @return A \code{draws} object from the
#'   \pkg{\link[posterior:posterior-package]{posterior}} package. See the
#'   \pkg{posterior} package documentation and vignettes for details on working
#'   with these objects.
#' 
#' @examples
#' fit <- stan_glm(mpg ~ wt + as.factor(cyl), data = mtcars)
#' as_draws_matrix(fit) # matrix format combines all chains 
#' as_draws_df(fit, regex_pars = "cyl")
#' posterior::summarize_draws(as_draws_array(fit))
#'
NULL

#' @rdname stanreg-draws-formats
#' @importFrom posterior as_draws
#' @method as_draws stanreg
#' @export
#' @export as_draws
as_draws.stanreg <- function(x, ...) {
  as_draws_df(x, ...)
}

#' @rdname stanreg-draws-formats
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix stanreg
#' @export
#' @export as_draws_matrix
as_draws_matrix.stanreg <- function(x, ...) {
  posterior::as_draws_matrix(
    as.matrix.stanreg(x, ...)
  )
}

#' @rdname stanreg-draws-formats
#' @importFrom posterior as_draws_array
#' @method as_draws_array stanreg
#' @export
#' @export as_draws_array
as_draws_array.stanreg <- function(x, ...) {
  if (used.sampling(x)) {
    posterior::as_draws_array(
      as.array.stanreg(x, ...)
    )
  } else {
    stop("For models not fit using MCMC use 'as_draws_matrix' instead of 'as_draws_array'",
         call. = FALSE)
  }
}

#' @rdname stanreg-draws-formats
#' @importFrom posterior as_draws_df
#' @method as_draws_df stanreg
#' @export
#' @export as_draws_df
as_draws_df.stanreg <- function(x, ...) {
  posterior::as_draws_df(
    if (used.sampling(x)) {
      as.array.stanreg(x, ...)
    } else {
      as.matrix.stanreg(x, ...)
    }
  )
}

#' @rdname stanreg-draws-formats
#' @importFrom posterior as_draws_list
#' @method as_draws_list stanreg
#' @export
#' @export as_draws_list
as_draws_list.stanreg <- function(x, ...) {
  posterior::as_draws_list(
    if (used.sampling(x)) {
      as.array.stanreg(x, ...)
    } else {
      as.matrix.stanreg(x, ...)
    }
  )
}

#' @rdname stanreg-draws-formats
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars stanreg
#' @export
#' @export as_draws_rvars
as_draws_rvars.stanreg <- function(x, ...) {
  posterior::as_draws_rvars(
    if (used.sampling(x)) {
      as.array.stanreg(x, ...)
    } else {
      as.matrix.stanreg(x, ...)
    }
  )
}
