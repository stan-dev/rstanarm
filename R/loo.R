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

#' Leave-one-out (LOO) and K-fold cross-validation
#' 
#' For models fit using MCMC, compute approximate leave-one-out cross-validation
#' (LOO) or, less preferably, the Widely Applicable Information Criterion (WAIC)
#' using the \pkg{\link[=loo-package]{loo}} package. Exact K-fold
#' cross-validation is also available. Compare two or more models using the
#' \code{\link[loo]{compare}} function.
#' 
#' @aliases loo waic compare
#'
#' @export
#' @templateVar stanregArg x
#' @template args-stanreg-object
#' @template reference-loo
#' 
#' @inheritParams loo::loo
#' 
#' @return An object of class 'loo'. See the 'Value' section in 
#'   \code{\link[loo]{loo}} and \code{\link[loo]{waic}} for details on the
#'   structure of these objects. The object returned by \code{kfold} also 
#'   has class 'kfold' in addition to 'loo'.
#'   
#' @details
#' \subsection{\code{loo}}{
#' The \code{loo} method for stanreg objects provides an interface to
#' the \pkg{\link[=loo-package]{loo}} package for approximate leave-one-out 
#' cross-validation (LOO). The LOO Information Criterion (LOOIC) has the same 
#' purpose as the Akaike Information Criterion (AIC) that is used by 
#' frequentists. Both are intended to estimate the expected log predictive 
#' density (ELPD) for a new dataset. However, the AIC ignores priors and assumes
#' that the posterior distribution is multivariate normal, whereas the functions
#' from the \pkg{loo} package do not make this distributional assumption and 
#' integrate over uncertainty in the parameters. This only assumes that any one 
#' observation can be omitted without having a major effect on the posterior 
#' distribution, which can be judged using the diagnostic plot provided by the 
#' \code{\link[loo]{plot.loo}} method and the warnings provided by the 
#' \code{\link[loo]{print.loo}} method (see the \emph{How to Use the rstanarm 
#' Package} vignette for an example of this process).
#' }
#' \subsection{\code{reloo}}{
#' The \code{reloo} function is provided as a possible remedy when the 
#' diagnositcs from \pkg{loo} reveal problems stemming from the posterior's 
#' sensitivity to particular observations. Warnings about Pareto \eqn{k}
#' estimates indicate observations for which the approximation to LOO is 
#' problematic. The \code{threshold} argument to \code{reloo} can be used to set
#' the \eqn{k} value above which an observation is flagged. If there are \eqn{J}
#' observations with \eqn{k} estimates above \code{threshold} then the
#' \code{reloo} function will refit the original model \eqn{J} times, each time
#' leaving out one of the \eqn{J} problematic observations. The pointwise
#' contributions of these observations to the total ELPD can then be computed
#' directly and be subsituted for the previous estimates from these \eqn{J}
#' observations that are stored in \code{loo_x} object.
#' }
#' \subsection{kfold}{
#' The \code{kfold} function performs exact K-fold cross-validation. The
#' \code{compare} function is also compatible with the objects returned by
#' \code{kfold}.
#' }
#'   
#' @seealso 
#' \code{\link[loo]{compare}} for comparing two or more models on LOO, WAIC, or
#' K-fold CV.
#' 
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
#' 
#' # compare on LOOIC
#' (loo1 <- loo(fit1))
#' loo2 <- loo(fit2)
#' compare(loo1, loo2)
#' plot(loo2)
#' 
#' # 10-fold cross-validation
#' (kfold1 <- kfold(fit1, K = 10))
#' kfold2 <- kfold(fit2, K = 10)
#' compare(kfold1, kfold2)
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
#' # compare on LOOIC
#' (loo_lalonde1 <- loo(lalonde1))
#' (loo_lalonde2 <- loo(lalonde2))
#' plot(loo_lalonde2, label_points = TRUE)
#' compare(loo_lalonde1, loo_lalonde2)
#' }
#' 
#' @importFrom loo loo loo.function compare
#' 
loo.stanreg <- function(x, ...) {
  if (!used.sampling(x)) 
    STOP_sampling_only("loo")
  loo.function(ll_fun(x), args = ll_args(x), ...)
}


#' @rdname loo.stanreg
#' @export
#' @param K The number of subsets of equal (if possible) size into which the 
#'   data will be randomly partitioned for performing K-fold cross-validation.
#'   The model is refit \code{K} times, each time leaving out one of the
#'   \code{K} subsets. If \code{K} is equal to the total number of observations
#'   in the data then K-fold cross-validation is equivalent to exact
#'   leave-one-out cross-validation.
#'
kfold <- function(x, K) {
  validate_stanreg_object(x)
  if (!used.sampling(x)) 
    STOP_sampling_only("kfold")
  stopifnot(!is.null(x$data), nrow(x$data) >= K)
  
  d <- x$data
  N <- nrow(d)
  perm <- sample.int(N)
  idx <- ceiling(seq(from = 1, to = N, length.out = K + 1))
  bin <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
  
  lppds <- list()
  for (k in 1:K) {
    message("Fitting model ", k, " out of ", K)
    omitted <- which(bin == k)
    fit_k <- update(x, data = d[-omitted, ], refresh = 0)
    lppds[[k]] <- log_lik(fit_k, newdata = d[omitted, ])
  }
  elpds <- unlist(lapply(lppds, function(x) {
    apply(x, 2, log_mean_exp)
  }))
  
  out <- list(
    elpd_kfold = sum(elpds),
    se_elpd_kfold = sqrt(N * var(elpds)),
    pointwise = cbind(elpd_kfold = elpds)
  )
  structure(out, class = c("kfold", "loo"), K = K)
}

print.kfold <- function(x, digits = 1, ...) {
  cat("\n", paste0(attr(x, "K"), "-fold"), "cross-validation\n\n")
  out <- data.frame(Estimate = x$elpd_kfold, SE = x$se_elpd_kfold, 
                    row.names = "elpd_kfold")
  .printfr(out, digits)
  invisible(x)
}
log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
}


#' @rdname loo.stanreg
#' @export 
#' @param loo_x The object returned by \code{loo} when applied to the stanreg 
#'   object \code{x}. If \code{loo_x} is missing then \code{reloo} runs 
#'   \code{loo} internally to create \code{loo_x} before proceeding. See Details
#'   for more information.
#' @param threshold Threshold for flagging estimates of the Pareto shape 
#'   parameters \eqn{k} estimated by \code{loo}. See Details.
#' 
reloo <- function(x, loo_x, threshold = 0.5, ...) {
  validate_stanreg_object(x)
  stopifnot(!is.null(x$data))
  
  if (missing(loo_x)) {
    loo_x <- loo(x, ...)
  } else {
    stopifnot(inherits(loo_x, "loo"))
    if (is.null(loo_x$pareto_k))
      stop("No Pareto k estimates found in 'loo' object.")
  }
  
  obs <- loo::pareto_k_ids(loo_x, threshold = threshold)
  if (!length(obs)) {
    message("No problematic observations found. Returning loo object.")
    return(loo_x)
  } 
  
  J <- length(obs)
  d <- x$data
  lls <- vector("list", J)
  
  message(
    J, " problematic observation(s) found.", 
    "\nModel will be refit ", J, " times."
  )
  for (j in 1:J) {
    message(
      "\nFitting model ", j, " out of ", J,
      " (leaving out observation ", obs[j], ")"
    )
    omitted <- obs[j]
    fit_j <- update(x, data = d[-omitted, ], refresh = 0)
    lls[[j]] <- log_lik(fit_j, newdata = d[omitted, ])
  }
  
  # replace parts of loo_x
  sel <- c("elpd_loo", "looic")
  elppds <- unlist(lapply(lls, log_mean_exp))
  loo_x$pointwise[obs, sel] <- cbind(elppds, -2 * elppds)
  loo_x[sel] <- with(loo_x, colSums(pointwise[, sel]))
  loo_x[paste0("se_", sel)] <- with(loo_x, {
    N <- nrow(pointwise)
    sqrt(N * apply(pointwise[, sel], 2, var))
  })
  
  # what should we do about pareto k's? for now setting them to 0
  loo_x$pareto_k[obs] <- 0
  loo_x
}


#' @rdname loo.stanreg
#' @export
#' @importFrom loo waic waic.function
#' @note The \code{...} is ignored for \code{waic}.
#' 
waic.stanreg <- function(x, ...) {
  if (!used.sampling(x)) 
    STOP_sampling_only("waic")
  waic.function(ll_fun(x), args = ll_args(x))
}



# function for loo.function -----------------------------------------------
# returns log-likelihood function for loo.function() and waic.function()
ll_fun <- function(x) {
  validate_stanreg_object(x)
  f <- family(x)
  if (!is(f, "family") || is_scobit(x))
    return(.ll_polr_i)
  
  get(paste0(".ll_", f$family, "_i"))
}


# arguments for loo.function ----------------------------------------------
# returns 'args' argument for loo.function() and waic.function()
ll_args <- function(object, newdata = NULL) {
  validate_stanreg_object(object)
  f <- family(object)
  draws <- nlist(f)
  has_newdata <- !is.null(newdata)
  if (has_newdata) {
    ppdat <- pp_data(object, as.data.frame(newdata))
    tmp <- pp_eta(object, ppdat)
    eta <- tmp$eta
    stanmat <- tmp$stanmat
    x <- ppdat$x
    y <- eval(formula(object)[[2L]], newdata)
  } else {
    stanmat <- as.matrix.stanreg(object)
    x <- get_x(object)
    y <- get_y(object)
  }
  
  if (is(f, "family") && !is_scobit(object)) {
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
    
  } else {
    stopifnot(is(object, "polr"))
    y <- as.integer(y)
    if (has_newdata) 
      x <- .validate_polr_x(object, x)
    data <- data.frame(y, x)
    draws$beta <- stanmat[, colnames(x), drop = FALSE]
    patt <- if (length(unique(y)) == 2L) "(Intercept)" else "|"
    zetas <- grep(patt, colnames(stanmat), fixed = TRUE, value = TRUE)
    draws$zeta <- stanmat[, zetas, drop = FALSE]
    draws$max_y <- max(y)
    if ("alpha" %in% colnames(stanmat)) {
      draws$alpha <- stanmat[, "alpha"]
      draws$f <- object$method
    }
  }
  
  data$offset <- object$offset
  if (!all(object$weights == 1)) 
    data$weights <- object$weights
  
  if (is.mer(object)) {
    b <- stanmat[, b_names(colnames(stanmat)), drop = FALSE]
    if (has_newdata) {
      Z_names <- ppdat$Z_names
      if (is.null(Z_names)) {
        b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
      } else {
        b <- pp_b_ord(b, Z_names)
      }
      z <- t(ppdat$Zt)
    } else {
      z <- get_z(object)
    }
    data <- cbind(data, as.matrix(z))
    draws$beta <- cbind(draws$beta, b)
  }
  
  nlist(data, draws, S = NROW(draws$beta), N = nrow(data))
}


# check intercept for polr models -----------------------------------------
# Check if a model fit with stan_polr has an intercept (i.e. if it's actually a 
# bernoulli model). If it doesn't have an intercept then the intercept column in
# x is dropped. This is only necessary if newdata is specified because otherwise
# the correct x is taken from the fitted model object.
.validate_polr_x <- function(object, x) {
  x0 <- get_x(object)
  has_intercept <- colnames(x0)[1L] == "(Intercept)" 
  if (!has_intercept && colnames(x)[1L] == "(Intercept)")
    x <- x[, -1L, drop = FALSE]
  x
}


# log-likelihood function helpers -----------------------------------------
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


# log-likelihood functions ------------------------------------------------
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
    if (y_i == 0) {
      val <- draws$alpha * log(linkinv(draws$zeta[, 1] - eta))
    } else if (y_i == 1) {
      val <- log1p(-linkinv(draws$zeta[, 1] - eta) ^ draws$alpha)
    } else {
      stop("Exponentiation only possible when there are exactly 2 outcomes.")
    }
  }
  .weighted(val, data$weights)
}
