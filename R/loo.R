# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

#' Leave-one-out cross-validation (LOO)
#' 
#' For models fit using MCMC, compute approximate leave-one-out cross-validation
#' (LOO) or, less preferably, the Widely Applicable Information Criterion (WAIC) 
#' using the \pkg{\link[=loo-package]{loo}} package. Compare two or more models 
#' using the \code{compare_models} function.
#' 
#' @aliases loo waic
#'
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template reference-loo
#' @inheritParams loo::loo
#' @return The \code{loo} and \code{waic} methods return an object of class
#'   'loo'. See the 'Value' sections in \code{\link[loo]{loo}} and
#'   \code{\link[loo]{waic}} for details on the structure of these objects.
#'   
#' @details 
#' The LOO Information Criterion (LOOIC) has the same purpose as the Aikaike
#' Information Criterion (AIC) that is used by frequentists. Both are intended
#' to estimate the expected log predicted density (ELPD) for a new dataset.
#' However, the AIC ignores priors and assumes that the posterior distribution
#' is multivariate normal, whereas the functions from the 
#' \pkg{\link[=loo-package]{loo}} package do not make this distributional 
#' assumption and integrate over uncertainty in the parameters. This only 
#' assumes that any one observation can be omitted without having a major effect
#' on the posterior distribution, which can be judged using the diagnostic plot
#' provided by the \code{\link[loo]{plot.loo}} method. The \emph{How to Use the
#' rstanarm Package} vignette has an example of this entire process.
#'   
#' @seealso 
#' \code{\link[loo]{loo-package}} (in particular the \emph{PSIS-LOO} section) 
#' for details on the computations implemented by the \pkg{loo} package and the 
#' interpretation of the Pareto \eqn{k} estimates displayed when using the 
#' \code{\link{plot.loo}} method.
#'   
#' @examples 
#' \dontrun{
#' SEED <- 42024
#' set.seed(SEED)
#' 
#' fit1 <- stan_glm(mpg ~ wt, data = mtcars, seed = SEED)
#' fit2 <- update(fit1, formula = . ~ . + cyl)
#' (loo1 <- loo(fit1))
#' loo2 <- loo(fit2)
#' compare_models(loos = list(loo1, loo2))
#' 
#' fit3 <- update(fit2, formula = . ~ . + am)
#' (comp <- compare_models(loos = list(loo1, loo2, loo(fit3))))
#' print(comp, digits = 2)
#' 
#' 
#' # Description of lalonde data at help("lalonde", package = "arm")
#' data(lalonde, package = "arm")
#' t7 <- student_t(df = 7) # prior for coefficients
#' 
#' f1 <- treat ~ re74 + re75 + educ + black + hisp + married + 
#'    nodegr + u74 + u75
#' lalonde1 <- stan_glm(f1, data = lalonde, family = binomial(link="logit"), 
#'                      prior = t7, cores = 4, seed = SEED)
#'                  
#' f2 <- treat ~ age + I(age^2) + educ + I(educ^2) + black + hisp + 
#'    married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) + u74 + u75   
#' lalonde2 <- update(lalonde1, formula = f2)
#' 
#' loos <- list(loo(lalonde1), loo(lalonde2))
#' compare_models(loos)
#' plot(loos[[2]], label_points = TRUE) # see help("plot.loo", package = "loo")
#' }
#' 
loo.stanreg <- function(x, ...) {
  if (!used.sampling(x)) 
    STOP_sampling_only("loo")
  if (!requireNamespace("digest", quietly = TRUE)) 
    stop("Please install the 'digest' package.")
  out <- loo.function(ll_fun(x$family), args = ll_args(x), ...)
  structure(out, family = family(x), name = deparse(substitute(x)),
            yhash = digest::digest(get_y(x), algo = "md5"))
}

#' @rdname loo.stanreg
#' @export
#' @note The \code{...} is ignored for \code{waic}.
#' 
waic.stanreg <- function(x, ...) {
  if (!used.sampling(x)) 
    STOP_sampling_only("waic")
  if (!requireNamespace("digest", quietly = TRUE)) 
    stop("Please install the 'digest' package.")
  out <- waic.function(ll_fun(x$family), args = ll_args(x))
  structure(out, family = family(x), name = deparse(substitute(x)), 
            yhash = digest::digest(get_y(x), algo = "md5"))
}

# returns log-likelihood function for loo() and waic()
ll_fun <- function(f) {
  if (is(f, "family")) {
    return(get(paste0(".ll_", f$family, "_i")))
  } else if (is.character(f)) {
    return(.ll_polr_i)
  } else {
    stop("'family' must be a family or a character string.", 
         call. = FALSE)
  }
}

# returns args argument for loo.function() and waic.function()
ll_args <- function(object) {
  f <- object$family
  draws <- nlist(f)
  stanmat <- as.matrix.stanreg(object)
  x <- get_x(object)
  y <- get_y(object)

  if (is(f, "family")) {
    fname <- f$family
    if (!is.binomial(fname)) {
      data <- data.frame(y, x)
    } else {
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L]
      } else {
        trials <- 1
        if (is.factor(y)) 
          y <- fac2bin(y)
        stopifnot(all(y %in% c(0, 1)))
      }
      data <- data.frame(y, trials, x)
    }
    draws$beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
    if (is.gaussian(fname)) 
      draws$sigma <- stanmat[, "sigma"]
    if (is.gamma(fname)) 
      draws$shape <- stanmat[, "shape"]
    if (is.ig(fname)) 
      draws$lambda <- stanmat[, "lambda"]
    if (is.nb(fname)) 
      draws$size <- stanmat[,"overdispersion"]
    
  } else if (is.character(f)) {
    stopifnot(is(object, "polr"))
    y <- as.integer(y)
    data <- data.frame(y, x)
    draws$beta <- stanmat[, colnames(x), drop = FALSE]
    zetas <- grep("|", colnames(stanmat), fixed = TRUE, value = TRUE)
    draws$zeta <- stanmat[, zetas, drop = FALSE]
    draws$max_y <- max(y)
    if ("alpha" %in% colnames(stanmat)) 
      draws$alpha <- stanmat[, "alpha"]
    
  } else {
    stop("'family' must be a family or a character string.", call. = FALSE)
  }
  
  data$offset <- object$offset
  if (!all(object$weights == 1)) 
    data$weights <- object$weights
  
  if (is.mer(object)) {
    z <- get_z(object)
    b <- stanmat[, b_names(colnames(stanmat)), drop = FALSE]
    data <- cbind(data, z)
    draws$beta <- cbind(draws$beta, b)
  }
  
  nlist(data, draws, S = NROW(draws$beta), N = nrow(data))
}


.xdata <- function(data) {
  sel <- c("y", "weights","offset", "trials")
  data[, -which(colnames(data) %in% sel)]
}
.mu <- function(data, draws) {
  eta <- as.vector(linear_predictor(draws$beta, .xdata(data), data$offset))
  draws$f$linkinv(eta)
}
.weighted <- function(val, w) {
  if (is.null(w)) {
    val
  } else {
    val * w
  } 
}

# log-likelihood functions
.ll_gaussian_i <- function(i, data, draws) {
  val <- dnorm(data$y, mean = .mu(data,draws), sd = draws$sigma, log = TRUE)
  .weighted(val, data$weights)
}
.ll_binomial_i <- function(i, data, draws) {
  val <- dbinom(data$y, size = data$trials, prob = .mu(data,draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_poisson_i <- function(i, data, draws) {
  val <- dpois(data$y, lambda = .mu(data, draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_neg_binomial_2_i <- function(i, data, draws) {
  val <- dnbinom(data$y, size = draws$size, mu = .mu(data, draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_Gamma_i <- function(i, data, draws) {
  val <- dgamma(data$y, shape = draws$shape, 
                rate = draws$shape / .mu(data,draws), log = TRUE)
  .weighted(val, data$weights)
}
.ll_inverse.gaussian_i <- function(i, data, draws) {
  mu <- .mu(data, draws)
  val <- 0.5 * log(draws$lambda / (2 * pi)) - 
         1.5 * log(data$y) -
         0.5 * draws$lambda * (data$y - mu)^2 / 
                        (data$y * mu^2)
  .weighted(val, data$weights)
}
.ll_polr_i <- function(i, data, draws) {
  eta <- linear_predictor(draws$beta, .xdata(data), data$offset)
  f <- draws$f
  J <- draws$max_y
  y_i <- data$y
  linkinv <- polr_linkinv(f)
  if (is.null(draws$alpha)) {
      if (y_i == 1) {
        val <- log(linkinv(draws$zeta[, 1] - eta))
      } else if (y_i == J) {
        val <- log1p(-linkinv(draws$zeta[, J-1] - eta))
      } else {
        val <- log(linkinv(draws$zeta[, y_i] - eta) - 
                     linkinv(draws$zeta[, y_i - 1L] - eta))
      }
  } else {
      if (y_i == 1) {
        val <- draws$alpha * log(linkinv(draws$zeta[, 1] - eta))
      } else if (y_i == J) {
        val <- log1p(-linkinv(draws$zeta[, J-1] - eta) ^ draws$alpha)
      } else {
        stop("Exponentiation only possible when there are exactly 2 outcomes.")
      }
  }
  .weighted(val, data$weights)
}


# Compare models
#
#' @rdname loo.stanreg
#' @export
#' @param loos A list of two or more objects of class "loo" returned by the
#'   \code{\link[=loo.stanreg]{loo}} method for 
#'   \code{\link[=stanreg-objects]{stanreg}} objects. See Examples.
#'   
#' @details 
#' \code{compare_models} is a wrapper around \code{\link[loo]{compare}} 
#' (\pkg{loo}) that performs some extra checks to make sure the models are
#' suitable for comparison.
#' 
#' @return \code{compare_models} returns a vector or matrix with class
#'   'compare.loo'. If \code{loos} contains more than two objects then a matrix
#'   is returned. This matrix summarizes the objects and also reports model
#'   weights (the posterior probability that each model has the best expected 
#'   out-of-sample predictive accuracy). If \code{loos} contains exactly two 
#'   objects then the difference in expected predictive accuracy and the 
#'   standard error of the difference are returned in addition to model weights.
#'   See the Details section in \code{\link[loo]{compare}}.
#' 
compare_models <- function(loos) {
  loos <- validate_loos(loos)
  comp <- do.call(compare, loos)
  if (!is.matrix(comp))  # will happen if there are only two models
    return(comp)
  
  stats <- c("looic", "se_looic", "elpd_loo", "se_elpd_loo", 
             "p_loo", "se_p_loo", "weights")
  structure(comp, dimnames = list(names(loos), stats))
}

validate_loos <- function(loos = list()) {
  if (!is.list(loos))
    stop("'loos' should be a list.", call. = FALSE)
  if (length(loos) <= 1)
    stop("'loos' should contain at least two objects.", call. = FALSE)
  families <- lapply(loos, attr, which = "family")
  ys <- lapply(loos, attr, which = "yhash")
  family_check <- sapply(families, function(x) {
    isTRUE(all.equal(x, families[[1]]))
  })
  ys_check <- sapply(ys, function(x) {
    isTRUE(all.equal(x, ys[[1]]))
  })
  if (!all(family_check))
    stop("Not all models have the same family/link.", call. = FALSE)
  if (!all(ys_check))
    stop("Not all models fit to the same y variable", call. = FALSE)
  setNames(loos, nm = lapply(loos, attr, which = "name"))
}
