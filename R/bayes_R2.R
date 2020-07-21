#' Compute a Bayesian version of R-squared or LOO-adjusted R-squared for
#' regression models.
#'
#' @aliases bayes_R2
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param re.form For models with group-level terms, \code{re.form} is
#'   passed to \code{\link{posterior_epred}} if specified.
#' @param ... Currently ignored.
#' 
#' @return A vector of R-squared values with length equal to the posterior
#'   sample size (the posterior distribution of R-squared).
#'   
#' @references
#' Andrew Gelman, Ben Goodrich, Jonah Gabry, and Aki Vehtari (2018). R-squared
#' for Bayesian regression models. \emph{The American Statistician}, to appear.
#' DOI: 10.1080/00031305.2018.1549100.
#' (\href{https://doi.org/10.1080/00031305.2018.1549100}{Journal},
#' \href{http://www.stat.columbia.edu/~gelman/research/published/bayes_R2_v3.pdf}{Preprint},
#' \href{https://avehtari.github.io/bayes_R2/bayes_R2.html}{Notebook})
#' 
#' @examples
#' fit <- stan_glm(
#'   mpg ~ wt + cyl, 
#'   data = mtcars, 
#'   QR = TRUE, 
#'   chains = 2, 
#'   refresh = 0
#' )
#' rsq <- bayes_R2(fit)
#' print(median(rsq))
#' hist(rsq)
#' 
#' loo_rsq <- loo_R2(fit)
#' print(median(loo_rsq))
#' 
#' # multilevel binomial model
#' if (!exists("example_model")) example(example_model)
#' print(example_model)
#' median(bayes_R2(example_model))
#' median(bayes_R2(example_model, re.form = NA)) # exclude group-level
#' 
bayes_R2.stanreg <- function(object, ..., re.form = NULL) {
    
    if (!used.sampling(object))
      STOP_sampling_only("bayes_R2")
    if (is_polr(object))
      stop("bayes_R2 is not available for stan_polr models.")
    
    fam <- family(object)$family
    if (!fam %in% c("gaussian", "binomial")) {
      stop("bayes_R2 is only available for Gaussian and binomial models.")
    }
    
    mu_pred <- posterior_epred(object, re.form = re.form)
    if (is.binomial(fam)) {
      y <- get_y(object)
      if (NCOL(y) == 2) {
        trials <- rowSums(y)
        trials_mat <- matrix(trials, nrow = nrow(mu_pred), ncol = ncol(mu_pred), 
                             byrow = TRUE)
        tmp <- mu_pred * trials_mat
        sigma2 <- rowMeans(tmp * (1 - mu_pred))
        mu_pred <- tmp
      } else {
        sigma2 <- rowMeans(mu_pred * (1 - mu_pred))
      }
    } else {
      sigma2 <- drop(as.matrix(object, pars = "sigma"))^2
    }
    
    var_mu_pred <- apply(mu_pred, 1, var)
    r_squared <- var_mu_pred / (var_mu_pred + sigma2)
    return(r_squared)
  }


#' @rdname bayes_R2.stanreg
#' @aliases loo_R2
#' @importFrom rstantools loo_R2
#' @export
#' 
loo_R2.stanreg <- function(object, ...) {
  if (!used.sampling(object))
    STOP_sampling_only("loo_R2")
  if (is_polr(object))
    stop("loo_R2 is not available for stan_polr models.")
  
  fam <- family(object)$family
  if (!fam %in% c("gaussian", "binomial")) {
    stop("loo_R2 is only available for Gaussian and binomial models.")
  }
  
  y <- get_y(object)
  log_ratios <- -log_lik(object)
  psis_object <- object[["loo"]][["psis_object"]]
  if (is.null(psis_object)) {
    psis_object <- loo::psis(log_ratios, r_eff = NA)
  }
  
  mu_pred <- posterior_epred(object)
  if (is.binomial(fam)) {
    if (is.factor(y)) {
      y <- fac2bin(y)
    } else if (NCOL(y) == 2) {
      trials <- rowSums(y)
      y <- y[, 1]
      trials_mat <- matrix(trials, nrow = nrow(mu_pred), ncol = ncol(mu_pred), 
                           byrow = TRUE)
      mu_pred <- mu_pred * trials_mat
    }
  }
  mu_pred_loo <- loo::E_loo(mu_pred, psis_object, log_ratios = log_ratios)$value
  err_loo <- mu_pred_loo - y
  
  S <- nrow(mu_pred)
  N <- ncol(mu_pred)

  # set the random seed as the seed used in the first chain and ensure
  # the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(object$stanfit@stan_args[[1]]$seed)

  # dirichlet weights 
  exp_draws <- matrix(rexp(S * N, rate = 1), nrow = S, ncol = N)
  wts <- exp_draws / rowSums(exp_draws)
  
  var_y <- (rowSums(sweep(wts, 2, y^2, FUN = "*")) -
            rowSums(sweep(wts, 2, y, FUN = "*"))^2) * (N/(N-1))
  
  var_err_loo <- (rowSums(sweep(wts, 2, err_loo^2, FUN = "*")) -
                  rowSums(sweep(wts, 2, err_loo, FUN = "*")^2)) * (N/(N-1))
  
  loo_r_squared <- 1 - var_err_loo / var_y
  loo_r_squared[loo_r_squared < -1] <- -1
  loo_r_squared[loo_r_squared > 1] <- 1
  return(loo_r_squared)
}


# internal ----------------------------------------------------------------
get_y_new <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    get_y(object)
  } else {
    eval(formula(object)[[2]], newdata)
  }
}
