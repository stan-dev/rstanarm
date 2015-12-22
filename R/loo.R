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
#' For models fit using MCMC (\code{algorithm="sampling"}), compute approximate
#' leave-one-out cross-validation (LOO) or the Widely Applicable Information
#' Criterion (WAIC) using the \pkg{\link[=loo-package]{loo}} package. Compare
#' two or more models using the \code{\link[loo]{compare}} function.
#' 
#' @aliases loo waic compare
#'
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template reference-loo
#' @inheritParams loo::loo
#' @return An object of class 'loo'. See \code{\link[loo]{loo}} and
#'   \code{\link[loo]{waic}}.
#'   
#' @seealso \code{\link[loo]{loo-package}}, \code{\link[loo]{compare}}, 
#'   \code{\link[loo]{plot.loo}}
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
#' compare(loo1, loo2)
#' plot(loo2)
#' 
#' 
#' # dataset description at help("lalonde", package = "arm")
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
#' (loo_lalonde1 <- loo(lalonde1))
#' (loo_lalonde2 <- loo(lalonde2))
#' plot(loo_lalonde2, label_points = TRUE)
#' compare(loo_lalonde1, loo_lalonde2)
#' }
#' 
#' @importFrom loo loo loo.function compare
#' 
loo.stanreg <- function(x, ...) {
  if (!used.sampling(x)) STOP_sampling_only("loo")
  loo.function(.llfun(x$family), args = .llargs(x), ...)
}

#' @rdname loo.stanreg
#' @export
#' @importFrom loo waic waic.function
#' @note The \code{...} is ignored for \code{waic}.
#' 
waic.stanreg <- function(x, ...) {
  if (!used.sampling(x)) STOP_sampling_only("waic")
  waic.function(.llfun(x$family), args = .llargs(x))
}

# returns log-likelihood function for loo() and waic()
.llfun <- function(f) {
  if (is(f, "family")) get(paste0(".ll_", f$family, "_i"))
  else if (is.character(f)) .ll_polr_i
  else stop("'family' must be a family or a character string")
}

# returns args argument for loo.function() and waic.function()
.llargs <- function(object) {
  f <- object$family
  draws <- nlist(f)
  stanmat <- as.matrix.stanreg(object)
  x <- get_x(object)
  y <- get_y(object)

  if (is(f, "family")) {
    if (!is.binomial(f$family)) data <- data.frame(y, x)
    else {
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L]
      } else {
        trials <- 1
        if (is.factor(y)) 
          y <- y != levels(y)[1L]
        stopifnot(all(y %in% c(0, 1)))
      }
      data <- data.frame(y, trials, x)
    }
    draws$beta <- stanmat[, 1:ncol(x), drop = FALSE]
    if (is.gaussian(f$family)) draws$sigma <- stanmat[, "sigma"]
    if (is.gamma(f$family)) draws$shape <- stanmat[, "shape"]
    if (is.ig(f$family)) draws$lambda <- stanmat[, "lambda"]
    if (is.nb(f$family)) draws$size <- stanmat[,"overdispersion"]
  }
  else if (is.character(f)) {
    stopifnot(is(object, "polr"))
    y <- as.integer(y)
    data <- data.frame(y, x)
    draws$beta <- stanmat[, colnames(x), drop = FALSE]
    draws$zeta <- stanmat[, grep("|", colnames(stanmat), fixed = TRUE, value = TRUE),
                          drop = FALSE]
    draws$max_y <- max(y)
  }
  else stop("'family' must be a family or a character string")
  
  data$offset <- object$offset
  if (!all(object$weights == 1)) data$weights <- object$weights
  
  if (is(object, "lmerMod")) {
    z <- get_z(object)
    b <- stanmat[, .bnames(colnames(stanmat)), drop = FALSE]
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
  if (is.null(w)) val
  else val * w
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
  val <- 
    if (y_i == 1) log(linkinv(draws$zeta[,1] - eta))
    else if (y_i == J) log1p(-linkinv(draws$zeta[,J-1] - eta))
    else log(linkinv(draws$zeta[,y_i] - eta) - 
               linkinv(draws$zeta[,y_i - 1L] - eta))
  
  .weighted(val, data$weights)
}
