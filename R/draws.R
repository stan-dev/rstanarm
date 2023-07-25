#' Create a \code{draws} object from a \code{stanreg} object
#'
#' Convert a \code{stanreg} object to a format supported by the \pkg{posterior}
#' package. To subset iterations, chains, or draws, use
#' \code{\link[posterior:subset_draws]{subset_draws}} after making the
#' \code{draws} object.
#'
#' @name draws-stanreg
#' @aliases as_draws as_draws_matrix as_draws_array as_draws_df as_draws_rvars as_draws_list
#'
#' @param x A \code{stanreg} object.
#' @param pars,regex_pars See \code{\link{as.matrix.stanreg}}.
#' @param ... Arguments passed to individual methods. 
#'
#' @examples
#' fit <- stan_glm(mpg ~ wt + cyl, data = mtcars)
#' head(as_draws_df(fit))
#' posterior <- as_draws_array(fit)
#' posterior::summarize_draws(posterior)
#'
#' as_draws_array(fit, variable = c("wt", "sigma"))
#'
NULL

#' @rdname draws-stanreg
#' @importFrom posterior as_draws
#' @method as_draws stanreg
#' @export
#' @export as_draws
as_draws.stanreg <- function(x, pars = NULL, regex_pars = NULL, ...) {
  as_draws_array(x, pars = pars, regex_pars = regex_pars, ...)
}

#' @rdname draws-stanreg
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix stanreg
#' @export
#' @export as_draws_matrix
as_draws_matrix.stanreg <- function(x, pars = NULL, regex_pars = FALSE, ...) {
  posterior::as_draws_matrix(
    as.matrix.stanreg(x, pars = pars, regex_pars = regex_pars, ...)
  )
}

#' @rdname draws-stanreg
#' @importFrom posterior as_draws_array
#' @method as_draws_array stanreg
#' @export
#' @export as_draws_array
as_draws_array.stanreg <- function(x, pars = NULL, regex_pars = FALSE, ...) {
  posterior::as_draws_array(
    as.array.stanreg(x, pars = pars, regex_pars = regex_pars, ...)
  )
}

#' @rdname draws-stanreg
#' @importFrom posterior as_draws_df
#' @method as_draws_df stanreg
#' @export
#' @export as_draws_df
as_draws_df.stanreg <- function(x, pars = NULL, regex_pars = FALSE, ...) {
  posterior::as_draws_df(
    as.array.stanreg(x, pars = pars, regex_pars = regex_pars, ...)
  )
}

#' @rdname draws-stanreg
#' @importFrom posterior as_draws_list
#' @method as_draws_list stanreg
#' @export
#' @export as_draws_list
as_draws_list.stanreg <- function(x, pars = NULL, regex_pars = FALSE, ...) {
  posterior::as_draws_list(
    as.array.stanreg(x, pars = pars, regex_pars = regex_pars, ...)
  )
}

#' @rdname draws-stanreg
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars stanreg
#' @export
#' @export as_draws_rvars
as_draws_rvars.stanreg <- function(x, pars = NULL, regex_pars = FALSE, ...) {
  posterior::as_draws_rvars(
    as.array.stanreg(x, pars = pars, regex_pars = regex_pars, ...)
  )
}
