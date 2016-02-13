# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

#' Draw from posterior predictive distribution
#' 
#' The posterior predictive distribution is the distribution of the outcome 
#' implied by the model after using the observed data to update our beliefs 
#' about the unknown parameters in the model. Simulating data from the posterior
#' predictive distribution using the observed predictors is useful for checking 
#' the fit of the model. Drawing from the posterior predictive distribution at 
#' interesting values of the predictors also lets us visualize how a 
#' manipulation of a predictor affects (a function of) the outcome(s). With new 
#' observations of predictor variables we can use posterior predictive 
#' distribution to generate predicted outcomes.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used. If \code{newdata} 
#'   is provided and any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdata}. This only applies if variables were transformed before 
#'   passing the data to one of the modeling functions and \emph{not} if 
#'   transformations were specified inside the model formula. Also see the Note
#'   section below for a note about using the \code{newdata} argument with with
#'   binomial models.
#' @param draws An integer indicating the number of draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param re.form If \code{object} contains \code{\link[=stan_glmer]{group-level}}
#'   parameters, a formula indicating which group-level parameters to 
#'   condition on when making predictions. \code{re.form} is specified in the 
#'   same form as for \code{\link[lme4]{predict.merMod}}. The default, 
#'   \code{NULL}, indicates that all estimated group-level parameters are 
#'   conditioned on. To refrain from conditioning on any group-level parameters,
#'   specify \code{NA} or \code{~0}. The \code{newdata} argument may include new
#'   \emph{levels} of the grouping factors that were specified when the model 
#'   was estimated, in which case the resulting posterior predictions 
#'   marginalize over the relevant variables.
#' @param fun An optional function to apply to the results. \code{fun} is found 
#'   by a call to \code{\link{match.fun}} and so can be specified as a function
#'   object, a string naming a function, etc.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently unused.
#' 
#' @return A \code{draws} by \code{nrow(newdata)} matrix of simulations
#'   from the posterior predictive distribution. Each row of the matrix is a
#'   vector of predictions generated using a single draw of the model parameters
#'   from the posterior distribution.
#'   
#' @note For binomial models with a number of trials greater than one (i.e., not
#'   Bernoulli models), if \code{newdata} is specified then it must include all 
#'   variables needed for computing the number of binomial trials to use for the
#'   predictions. For example if the left-hand side of the model formula is 
#'   \code{cbind(successes, failures)} then both \code{successes} and 
#'   \code{failures} must be in \code{newdata}. The particular values of 
#'   \code{successes} and \code{failures} in \code{newdata} do not matter so 
#'   long as their sum is the desired number of trials. If the left-hand side of
#'   the model formula were \code{cbind(successes, trials - successes)} then
#'   both \code{trials} and \code{successes} would need to be in \code{newdata},
#'   probably with \code{successes} set to \code{0} and \code{trials} specifying
#'   the number of trials. See the Examples section below and the 
#'   \emph{How to Use the rstanarm Package} for examples.
#' 
#' @seealso \code{\link{pp_check}} for graphical posterior predictive checks.
#'   Examples of posterior predictive checking can also be found in the
#'   \pkg{rstanarm} vignettes and demos.
#'   
#' @examples
#' yrep <- posterior_predict(example_model)
#' table(yrep)
#' 
#' \dontrun{
#' # Using newdata
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' fit3 <- stan_glm(counts ~ outcome + treatment, family = poisson(link="log"),
#'                 prior = normal(0, 1), prior_intercept = normal(0, 5))
#' nd <- data.frame(treatment = factor(rep(1,3)), outcome = factor(1:3))
#' ytilde <- posterior_predict(fit3, nd, draws = 500)
#' print(dim(ytilde))  # 500 by 3 matrix (draws by nrow(nd))
#' ytilde <- data.frame(count = c(ytilde), 
#'                      outcome = rep(nd$outcome, each = 500))
#' ggplot(ytilde, aes(x=outcome, y=count)) + 
#'   geom_boxplot() + 
#'   ylab("predicted count")
#' 
#' 
#' # Using newdata with a binomial model
#' # example_model is binomial so we need to set
#' # the number of trials to use for prediction.
#' # This could be a different number for each 
#' # row of newdata or the same for all rows.
#' # Here we'll use the same value for all.
#' nd <- lme4::cbpp
#' print(formula(example_model))  # cbind(incidence, size - incidence) ~ ...
#' nd$size <- max(nd$size) + 1L   # number of trials
#' nd$incidence <- 0  # set to 0 so size - incidence = number of trials
#' ytilde <- posterior_predict(example_model, newdata = nd)
#' 
#' 
#' # Using fun argument to transform predictions
#' mtcars2 <- mtcars
#' mtcars2$log_mpg <- log(mtcars2$mpg)
#' fit <- stan_glm(log_mpg ~ wt, data = mtcars2)
#' ytilde <- posterior_predict(fit, fun = exp)
#' }
#' 
posterior_predict <- function(object, newdata = NULL, draws = NULL, 
                              re.form = NULL, fun = NULL, seed = NULL, ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object.")
  if (used.optimizing(object))
    STOP_not_optimizing("posterior_predict")
  if (!is.null(seed)) 
    set.seed(seed)
  if (!is.null(fun)) 
    fun <- match.fun(fun)
  if (!is.null(newdata)) {
    if ("gam" %in% names(object))
      stop("'posterior_predict' with 'newdata' not yet supported ", 
           "for models estimated via 'stan_gamm4'.")
    newdata <- as.data.frame(newdata)
    if (any(is.na(newdata))) 
      stop("Currently NAs are not allowed in 'newdata'.")
  }
  dat <- pp_data(object, newdata, re.form, ...)
  ppargs <- pp_args(object, data = pp_eta(object, dat, draws))
  if (!is(object, "polr") && is.binomial(family(object)$family))
    ppargs$trials <- pp_binomial_trials(object, newdata)
  
  ppfun <- pp_fun(object)
  ytilde <- do.call(ppfun, ppargs)
  if (!is.null(newdata) && nrow(newdata) == 1L) 
    ytilde <- t(ytilde)
  if (!is.null(fun)) 
    ytilde <- do.call(fun, list(ytilde))
  
  return(ytilde)
}


# functions to draw from the various posterior predictive distributions
pp_fun <- function(object) {
  suffix <- if (is(object, "polr")) "polr" else family(object)$family
  get(paste0(".pp_", suffix), mode = "function")
}

.pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s,], sigma[s])
  }))
}
.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s, ])
  }))
}
.pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s, ])
  }))
}
.pp_neg_binomial_2 <- function(mu, size) {
  t(sapply(1:nrow(mu), function(s) {
    rnbinom(ncol(mu), size = size[s], mu = mu[s, ])
  }))
}
.pp_Gamma <- function(mu, shape) {
  t(sapply(1:nrow(mu), function(s) {
    rgamma(ncol(mu), shape = shape[s], rate = shape[s] / mu[s, ])
  }))
}
.rinvGauss <- function(n, mu, lambda) {
  # draw from inverse gaussian distribution
  mu2 <- mu^2
  y <- rnorm(n)^2
  z <- runif(n)
  tmp <- (mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * y^2))
  x <- mu + tmp / (2 * lambda)
  ifelse(z <= (mu / (mu + x)), x, mu2 / x)
}
.pp_inverse.gaussian <- function(mu, lambda) {
  t(sapply(1:nrow(mu), function(s) {
    .rinvGauss(ncol(mu), mu = mu[s,], lambda = lambda[s])
  }))
}
.pp_polr <- function(eta, zeta, linkinv, alpha = NULL) {
  n <- ncol(eta)
  q <- ncol(zeta)
  if (!is.null(alpha)) {
    t(sapply(1:nrow(eta), FUN = function(s) {
      tmp <- matrix(zeta[s,], n, q, byrow = TRUE) - eta[s, ]
      pr <- matrix(linkinv(tmp)^alpha, , q)
      rbinom(ncol(eta), size = 1, prob = pr[s, ])
    }))
  } else {
    t(sapply(1:nrow(eta), FUN = function(s) {
      tmp <- matrix(zeta[s, ], n, q, byrow = TRUE) - eta[s, ]
      cumpr <- matrix(linkinv(tmp), , q)
      fitted <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
      apply(fitted, 1, function(p) which(rmultinom(1, 1, p) == 1))
    }))
  }
}


# create list of arguments to pass to the function returned by pp_fun
#
# @param object stanreg object
# @data output from pp_eta (named list with eta and stanmat)
# @return named list
pp_args <- function(object, data) {
  stanmat <- data$stanmat
  eta <- data$eta
  stopifnot(is.stanreg(object), is.matrix(stanmat))
  inverse_link <- linkinv(object)
  if (is(object, "polr")) {
    zeta <- stanmat[, grep("|", colnames(stanmat), value = TRUE, fixed = TRUE)]
    args <- nlist(eta, zeta, linkinv = inverse_link)
    if ("alpha" %in% colnames(stanmat))
      args$alpha <- stanmat[, "alpha"]
    return(args)
  }
  
  args <- list(mu = inverse_link(eta))
  famname <- family(object)$family
  if (is.gaussian(famname)) {
    args$sigma <- stanmat[, "sigma"]
  } else if (is.gamma(famname)) {
    args$shape <- stanmat[, "shape"]
  } else if (is.ig(famname)) {
    args$lambda <- stanmat[, "lambda"]
  } else if (is.nb(famname)) {
    args$size <- stanmat[, "overdispersion"]
  }
  args
}

# create eta and stanmat (matrix of posterior draws)
# 
# @param object stanreg object
# @param data output from pp_data()
# @param draws number of draws
# @return linear predictor "eta" and matrix of posterior draws stanmat"
pp_eta <- function(object, data, draws = NULL) {
  x <- data$x
  S <- posterior_sample_size(object)
  if (is.null(draws)) 
    draws <- S
  if (draws > S) {
    err <- paste0("'draws' should be <= posterior sample size (", 
                  S, ").")
    stop(err)
  }
  some_draws <- isTRUE(draws < S)
  if (some_draws)
    samp <- sample(S, draws)
  if (is.null(data$Zt)) {
    stanmat <- as.matrix.stanreg(object)
    beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
    if (some_draws) 
      beta <- beta[samp, , drop = FALSE]
    eta <- linear_predictor(beta, x, data$offset)
  } else {
    stanmat <- as.matrix(object$stanfit)
    beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
    if (some_draws) 
      beta <- beta[samp, , drop = FALSE]
    eta <- linear_predictor(beta, x, data$offset)
    b <- stanmat[, grepl("^b\\[", colnames(stanmat)), drop = FALSE]
    if (some_draws) 
      b <- b[samp, , drop = FALSE]
    if (is.null(data$Z_names)) {
      b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
    } else {
      b <- pp_b_ord(b, data$Z_names)
    }
    eta <- eta + as.matrix(b %*% data$Zt)
  }
  nlist(eta, stanmat)
}

pp_b_ord <- function(b, Z_names) {
  ord <- sapply(Z_names, FUN = function(x) {
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1) 
      return(m)
    if (len > 1) 
      stop("multiple matches bug")
    x <- sub(" (.*):.*$", " \\1:_NEW_\\1", x)
    grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
  })
  b[, ord, drop = FALSE]
}

# Number of trials for binomial models
pp_binomial_trials <- function(object, newdata) {
  y <- if (is.null(newdata))
    get_y(object) else eval(formula(object)[[2L]], newdata)
  if (NCOL(y) == 2L) 
    return(rowSums(y))
  rep(1, NROW(y))
}
