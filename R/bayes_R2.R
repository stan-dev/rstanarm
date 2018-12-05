#' Compute a Bayesian version of R-squared for regression models
#'
#' @aliases bayes_R2
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param newdata Similar to the \code{newdata} argument to 
#'   \code{\link{posterior_linpred}} and \code{\link{posterior_predict}} except,
#'   in addition to new observations of the predictors, new observations of the
#'   \emph{outcome} must be also included. See the \strong{Examples} section below.
#' @param re.form,offset For models with group-level terms, these arguments are 
#'   passed to \code{\link{posterior_linpred}} if the \code{newdata} argument is
#'   specified.
#' @param ... Currently ignored.
#' 
#' @return A vector of Bayesian R-squared values with length equal to the 
#'   posterior sample size.
#'   
#' @seealso \url{https://github.com/jgabry/bayes_R2}
#' 
#' @examples
#' fit <- stan_glm(mpg ~ wt + cyl, data = mtcars, QR = TRUE, chains = 2)
#' rsq <- bayes_R2(fit)
#' print(median(rsq))
#' 
#' # specifying newdata (including outcome variable 'mpg')
#' nd <- data.frame(mpg = c(10, 20, 30), wt = c(4, 3, 2), cyl = c(8, 6, 4))
#' rsq_new <- bayes_R2(fit, newdata = nd)
#' print(median(rsq_new))
#' 
#' # multilevel binomial model
#' if (!exists("example_model")) example(example_model)
#' print(example_model)
#' median(bayes_R2(example_model))
#' median(bayes_R2(example_model, re.form = NA)) # exclude group-level
#' 
bayes_R2.stanreg <-
  function(object,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           ...) {
    
    if (!used.sampling(object))
      STOP_sampling_only("bayes_R2")
    if (is_polr(object))
      stop("bayes_R2 is not available for stan_polr models.")
    
    fam <- family(object)$family
    if (!fam %in% c("gaussian", "binomial")) {
      stop("bayes_R2 is only available for Gaussian and binomial models.")
    }
    
    yhat <- posterior_linpred(
      object,
      transform = TRUE,
      newdata = newdata,
      re.form = re.form,
      offset = offset
    )
    
    if (is.binomial(fam)) {
      y <- get_y_new(object, newdata = newdata)
      if (NCOL(y) == 2) {
        trials <- rowSums(y)
        yhat <- yhat %*% diag(trials)
      }
      sigma2 <- rowMeans(yhat * (1 - yhat))
    } else {
      sigma2 <- drop(as.matrix(fit, pars = "sigma"))^2
    }
    
    var_yhat <- apply(yhat, 1, var)
    r_squared <- var_yhat / (var_yhat + sigma2)
    return(r_squared)
  }


# internal ----------------------------------------------------------------
# get_y_new <- function(object, newdata = NULL) {
#   if (is.null(newdata)) {
#     get_y(object)
#   } else {
#     eval(formula(object)[[2]], newdata)
#   }
# }
