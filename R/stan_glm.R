# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Bayesian generalized linear models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Generalized linear modeling with optional prior distributions for the
#' coefficients, intercept, and auxiliary parameters.
#'
#' @export
#' @templateVar armRef (Ch. 3-6)
#' @templateVar pkg stats
#' @templateVar pkgfun glm
#' @templateVar sameargs model,offset,weights 
#' @templateVar rareargs na.action,contrasts
#' @templateVar fun stan_glm, stan_glm.nb
#' @templateVar fitfun stan_glm.fit
#' @template return-stanreg-object
#' @template return-stanfit-object
#' @template see-also
#' @template args-formula-data-subset
#' @template args-same-as
#' @template args-same-as-rarely
#' @template args-dots
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_aux
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template args-sparse
#' @template reference-gelman-hill
#' @template reference-muth
#' 
#' @param family Same as \code{\link[stats]{glm}}, except negative binomial GLMs
#'   are also possible using the \code{\link{neg_binomial_2}} family object.
#' @param y In \code{stan_glm}, logical scalar indicating whether to
#'   return the response vector. In \code{stan_glm.fit}, a response vector.
#' @param x In \code{stan_glm}, logical scalar indicating whether to
#'   return the design matrix. In \code{stan_glm.fit}, usually a design matrix
#'   but can also be a list of design matrices with the same number of rows, in
#'   which case the first element of the list is interpreted as the primary design
#'   matrix and the remaining list elements collectively constitute a basis for a
#'   smooth nonlinear function of the predictors indicated by the \code{formula}
#'   argument to \code{\link{stan_gamm4}}.
#' @param mean_PPD A logical value indicating whether the sample mean of the
#'   posterior predictive distribution of the outcome should be calculated in
#'   the \code{generated quantities} block. If \code{TRUE} then \code{mean_PPD}
#'   is computed and displayed as a diagnostic in the
#'   \link[=print.stanreg]{printed output}. The default is \code{TRUE} except if
#'   \code{algorithm=="optimizing"}. A useful heuristic is to check if
#'   \code{mean_PPD} is plausible when compared to \code{mean(y)}. If it is
#'   plausible then this does \emph{not} mean that the model is good in general
#'   (only that it can reproduce the sample mean), but if \code{mean_PPD} is
#'   implausible then there may be something wrong, e.g., severe model
#'   misspecification, problems with the data and/or priors, computational
#'   issues, etc.
#' 
#' @details The \code{stan_glm} function is similar in syntax to 
#'   \code{\link[stats]{glm}} but rather than performing maximum likelihood 
#'   estimation of generalized linear models, full Bayesian estimation is 
#'   performed (if \code{algorithm} is \code{"sampling"}) via MCMC. The Bayesian
#'   model adds priors (independent by default) on the coefficients of the GLM.
#'   The \code{stan_glm} function calls the workhorse \code{stan_glm.fit}
#'   function, but it is also possible to call the latter directly.
#'   
#'   The \code{stan_glm.nb} function, which takes the extra argument 
#'   \code{link}, is a wrapper for \code{stan_glm} with \code{family = 
#'   \link{neg_binomial_2}(link)}.
#'   
#' @seealso The various vignettes for \code{stan_glm} at
#'   \url{http://mc-stan.org/rstanarm/articles/}.
#' 
#' @examples
#' if (!grepl("^sparc",  R.version$platform)) {
#' ### Linear regression
#' mtcars$mpg10 <- mtcars$mpg / 10
#' fit <- stan_glm(
#'   mpg10 ~ wt + cyl + am,            
#'   data = mtcars, 
#'   QR = TRUE,
#'   # for speed of example only (default is "sampling")
#'   algorithm = "fullrank",
#'   refresh = 0 
#'  ) 
#'                 
#' plot(fit, prob = 0.5)
#' plot(fit, prob = 0.5, pars = "beta")
#' plot(fit, "hist", pars = "sigma")
#' }
#' \donttest{
#' ### Logistic regression
#' head(wells)
#' wells$dist100 <- wells$dist / 100
#' fit2 <- stan_glm(
#'   switch ~ dist100 + arsenic, 
#'   data = wells, 
#'   family = binomial(link = "logit"), 
#'   prior_intercept = normal(0, 10),
#'   QR = TRUE,
#'   refresh = 0,
#'   # for speed of example only
#'   chains = 2, iter = 200 
#' )
#' print(fit2)
#' prior_summary(fit2)
#' 
#' # ?bayesplot::mcmc_areas
#' plot(fit2, plotfun = "areas", prob = 0.9,
#'      pars = c("(Intercept)", "arsenic"))
#' 
#' # ?bayesplot::ppc_error_binned
#' pp_check(fit2, plotfun = "error_binned") 
#' 
#' 
#' ### Poisson regression (example from help("glm")) 
#' count_data <- data.frame(
#'  counts = c(18,17,15,20,10,20,25,13,12),
#'  outcome = gl(3,1,9),
#'  treatment = gl(3,3)
#' )
#' fit3 <- stan_glm(
#'   counts ~ outcome + treatment, 
#'   data = count_data, 
#'   family = poisson(link="log"),
#'   prior = normal(0, 2),
#'   refresh = 0,
#'   # for speed of example only
#'   chains = 2, iter = 250 
#' ) 
#' print(fit3)
#' 
#' bayesplot::color_scheme_set("viridis")
#' plot(fit3)
#' plot(fit3, regex_pars = c("outcome", "treatment"))
#' plot(fit3, plotfun = "combo", regex_pars = "treatment") # ?bayesplot::mcmc_combo
#' posterior_vs_prior(fit3, regex_pars = c("outcome", "treatment"))
#' 
#' ### Gamma regression (example from help("glm"))
#' clotting <- data.frame(log_u = log(c(5,10,15,20,30,40,60,80,100)),
#'                        lot1 = c(118,58,42,35,27,25,21,19,18),
#'                        lot2 = c(69,35,26,21,18,16,13,12,12))
#' fit4 <- stan_glm(
#'   lot1 ~ log_u, 
#'   data = clotting, 
#'   family = Gamma(link="log"),
#'   iter = 500, # for speed of example only
#'   refresh = 0
#'  ) 
#' print(fit4, digits = 2)
#' 
#' fit5 <- update(fit4, formula = lot2 ~ log_u)
#' 
#' # ?bayesplot::ppc_dens_overlay
#' bayesplot::bayesplot_grid(
#'   pp_check(fit4, seed = 123), 
#'   pp_check(fit5, seed = 123),
#'   titles = c("lot1", "lot2")
#' ) 
#' 
#' 
#' ### Negative binomial regression
#' fit6 <- stan_glm.nb(
#'   Days ~ Sex/(Age + Eth*Lrn), 
#'   data = MASS::quine, 
#'   link = "log", 
#'   prior_aux = exponential(1.5, autoscale=TRUE),
#'   chains = 2, iter = 200, # for speed of example only
#'   refresh = 0
#' ) 
#' 
#' prior_summary(fit6)
#' bayesplot::color_scheme_set("brightblue")
#' plot(fit6)
#' pp_check(fit6, plotfun = "hist", nreps = 5) # ?bayesplot::ppc_hist
#' 
#' # 80% interval of estimated reciprocal_dispersion parameter
#' posterior_interval(fit6, pars = "reciprocal_dispersion", prob = 0.8)
#' plot(fit6, "areas", pars = "reciprocal_dispersion", prob = 0.8)
#' }
#'
stan_glm <-
  function(formula,
           family = gaussian(),
           data,
           weights,
           subset,
           na.action = NULL,
           offset = NULL,
           model = TRUE,
           x = FALSE,
           y = TRUE,
           contrasts = NULL,
           ...,
           prior = default_prior_coef(family),
           prior_intercept = default_prior_intercept(family),
           prior_aux = exponential(autoscale=TRUE),
           prior_PD = FALSE,
           algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
           mean_PPD = algorithm != "optimizing",
           adapt_delta = NULL,
           QR = FALSE,
           sparse = FALSE) {
    
  algorithm <- match.arg(algorithm)
  family <- validate_family(family)
  validate_glm_formula(formula)
  data <- validate_data(data, if_missing = environment(formula))
  
  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "subset", "weights", "na.action", "offset"), 
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mf <- check_constant_vars(mf)
  mt <- attr(mf, "terms")
  Y <- array1D_check(model.response(mf, type = "any"))
  if (is.empty.model(mt))
    stop("No intercept or predictors specified.", call. = FALSE)
  X <- model.matrix(mt, mf, contrasts)
  weights <- validate_weights(as.vector(model.weights(mf)))
  offset <- validate_offset(as.vector(model.offset(mf)), y = Y)
  if (binom_y_prop(Y, family, weights)) {
    y1 <- as.integer(as.vector(Y) * weights)
    Y <- cbind(y1, y0 = weights - y1)
    weights <- double(0)
  }
  
  if (prior_PD) {
    # can result in errors (e.g. from poisson) if draws from prior are weird
    mean_PPD <- FALSE
  }

  stanfit <- stan_glm.fit(
    x = X,
    y = Y,
    weights = weights,
    offset = offset,
    family = family,
    prior = prior,
    prior_intercept = prior_intercept,
    prior_aux = prior_aux,
    prior_PD = prior_PD,
    algorithm = algorithm,
    mean_PPD = mean_PPD,
    adapt_delta = adapt_delta,
    QR = QR,
    sparse = sparse,
    ...
  )
  if (algorithm != "optimizing" && !is(stanfit, "stanfit")) return(stanfit)
  if (family$family == "Beta regression") {
    family$family <- "beta"
  }

  sel <- apply(X, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  X <- X[ , !sel, drop = FALSE]  

  fit <- nlist(stanfit, algorithm, family, formula, data, offset, weights,
               x = X, y = Y, model = mf,  terms = mt, call, 
               na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"), 
               stan_function = "stan_glm")
  
  out <- stanreg(fit)
  if (algorithm == "optimizing") {
    out$log_p <- stanfit$log_p
    out$log_g <- stanfit$log_g
    out$psis <- stanfit$psis
    out$ir_idx <- stanfit$ir_idx
    out$diagnostics <- stanfit$diagnostics
  }
  out$compute_mean_PPD <- mean_PPD
  out$xlevels <- .getXlevels(mt, mf)
  if (!x) 
    out$x <- NULL
  if (!y) 
    out$y <- NULL
  if (!model) 
    out$model <- NULL
  
  return(out)
}

#' @rdname stan_glm
#' @export
#' @param link For \code{stan_glm.nb} only, the link function to use. See 
#'   \code{\link{neg_binomial_2}}.
#'   
stan_glm.nb <- 
  function(formula,
           data,
           weights,
           subset,
           na.action = NULL,
           offset = NULL,
           model = TRUE,
           x = FALSE,
           y = TRUE,
           contrasts = NULL,
           link = "log",
           ...,
           prior = default_prior_coef(family),
           prior_intercept = default_prior_intercept(family),
           prior_aux = exponential(autoscale=TRUE),
           prior_PD = FALSE,
           algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
           mean_PPD = algorithm != "optimizing",
           adapt_delta = NULL,
           QR = FALSE) {
    
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc[[1L]] <- quote(stan_glm)
  mc$link <- NULL
  mc$family <- neg_binomial_2(link = link)
  out <- eval(mc, parent.frame())
  out$call <- call
  out$stan_function <- "stan_glm.nb"
  return(out)
}
