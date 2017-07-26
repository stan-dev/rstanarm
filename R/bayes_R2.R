#' Compute a Bayesian version of R-squared for regression models
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param newdata Similar to the \code{newdata} argument to 
#'   \code{\link{posterior_linpred}} and \code{\link{posterior_predict}} except,
#'   in addition to new observations of the predictors, new observations of the
#'   \emph{outcome} must be also included.
#' @param re.form,offset For models with group-level terms, these arguments are 
#'   passed to \code{\link{posterior_linpred}} if the \code{newdata} argument is
#'   specified.
#' @param ... Currently ignored.
#' 
#' @return A vector of Bayesian R-squared values with length equal to the 
#'   posterior sample size.
#' 
#' @examples
#' fit <- stan_glm(mpg ~ wt + cyl, data = mtcars)
#' rsq <- bayes_R2(fit)
#' print(median(rsq))
#' 
#' nd <- data.frame(mpg = c(10, 20, 30), wt = c(4, 3, 2), cyl = c(8, 6, 4))
#' rsq_new <- bayes_R2(fit, newdata = nd)
#' print(median(rsq_new))
#' 
#' if (!exists("example_model")) example(example_model)
#' print(example_model)
#' median(bayes_R2(example_model))
#' median(bayes_R2(example_model, re.form = NA)) # exclude group-level
#' 
bayes_R2 <-
  function(object,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           ...) {
    
    if (!used.sampling(object))
      STOP_sampling_only("bayes_R2")
    if (is_polr(object))
      stop("Not available for stan_polr models.")
    
    y <- get_y_new(object, newdata = newdata)
    ypred <- posterior_linpred(
      object,
      transform = TRUE,
      newdata = newdata,
      re.form = re.form,
      offset = offset
    )
    
    if (is.binomial(family(object)$family)) {
      if (is.factor(y)) {
        y <- fac2bin(y)
      } else if (NCOL(y) == 2) {
        trials <- rowSums(y)
        y <- y[, 1]
        ypred <- ypred %*% diag(trials)
      }
    }
    
    e <- -1 * sweep(ypred, 2, y)
    var_ypred <- apply(ypred, 1, var)
    var_e <- apply(e, 1, var)
    var_ypred / (var_ypred + var_e)
  }


# internal ----------------------------------------------------------------
get_y_new <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    get_y(object)
  } else {
    eval(formula(object)[[2]], newdata)
  }
}
