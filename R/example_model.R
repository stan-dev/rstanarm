#' Example model
#' 
#' A pre-fit model for use in \pkg{rstanarm} examples. 
#' 
#' @name example_model
#' @format A \code{\link[=stanreg-objects]{stanreg}} object containing the
#'   output from fitting the model in the examples section.
#'   The chains and iter arguments are specified to make this example be 
#'   small in size. In practice, we recommend that they be left unspecified
#'   in order to use the default values (4 and 2000 respectively) or increased
#'   if there are convergence problems. The cores argument is optional and on
#'   a multicore system, the user may well want to set that equal to the number
#'   of chains being executed.
#' 
#' @seealso \code{\link[lme4]{cbpp}} for a description of the data.
#' @examples
#' \dontrun{ 
#' example_model <- 
#'   stan_glmer(cbind(incidence, size - incidence) ~ size + period + (1|herd),
#'              data = lme4::cbpp, family = binomial,
#'              chains = 2, cores = 1, seed = 12345, iter = 500)
#' }
NULL