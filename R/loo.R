#' Leave-one-out cross-validation (LOO)
#' 
#' Compute approximate leave-one-out cross-validation (LOO) or the Widely 
#' Applicable Information Criterion (WAIC) using the
#' \pkg{\link[=loo-package]{loo}} package.
#' 
#' @export
#' @inheritParams loo::loo
#' @param x A fitted model object returned by one of the \pkg{rstanarm} modeling
#'   functions. This will be a list with class 'stanreg' as well as at least one
#'   of 'lm', 'glm', 'polr', 'lmerMod', or 'aov'.
#' @return An object of class 'loo'. See \code{\link[loo]{loo}} and
#'   \code{\link[loo]{waic}}.
#' 
#' @seealso \code{\link[loo]{loo-package}}
#' @importFrom loo loo loo.function
#' 
#' 
loo.stanreg <- function(x, ...) {
  if (x$algorithm != "sampling")
    stop("Only available for MCMC (algorithm = 'sampling').", call. = FALSE)
  loo.function(.llfun(x), args = .llargs(x), ...)
}

#' @rdname loo.stanreg
#' @export
#' @importFrom loo waic waic.function
#' @note The \code{...} is ignored for \code{waic}.
#' 
waic.stanreg <- function(x, ...) {
  if (x$algorithm != "sampling")
    stop("Only available for MCMC (algorithm = 'sampling').", call. = FALSE)
  waic.function(.llfun(x), args = .llargs(x))
}

# returns args argument for loo() and waic()
.llargs <- function(x) {
  args <- list()
  args$family <- x$family
  stanmat <- as.matrix(x$stanfit)
  nms <- names(x)
  args$x <- if ("x" %in% nms) x$x else model.matrix(x)
  args$y <- if ("y" %in% nms) x$y else model.response(model.frame(x))
  args$offset <- x$offset
  args$weights <- if (!is.null(x$weights)) x$weights else 1
  if (is(args$family, "family")) {
    if (args$family$family == "gaussian") 
      args$sigma <- stanmat[, "sigma"]
    if (args$family$family == "Gamma")
      args$shape <- stanmat[, "shape"]
    if (args$family$family == "inverse.gaussian")
      args$lambda <- stanmat[, "lambda"]
  }
  if (is(x, "polr")) {
    args$beta <- stanmat[,colnames(args$x),drop = FALSE]
    args$zeta <- stanmat[,grep("|", colnames(stanmat), fixed = TRUE, value = TRUE),drop=FALSE]
  } else {
    args$beta <- stanmat[, 1:ncol(args$x)]
  }
  do.call(".llargs_prep", args)
}
.llargs_prep <- function(family, x, y, weights, beta,
                         # remaining arguments might be applicable
                         sigma = NULL, zeta = NULL, shape = NULL,
                         lambda = NULL, offset = NULL) {
  draws <- nlist(beta)
  f <- family
  draws$f <- f
  if (is(f, "family")) {
    if (f$family != "binomial") data <- data.frame(y, x)
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
  }
  else if (is.character(f)) {
    y <- as.integer(y)
    data <- data.frame(y,x)
    draws$max_y <- max(y)
  }
  else stop("'family' must be a family or a character string")
  
  if (!is.null(sigma))  draws$sigma <- sigma
  if (!is.null(zeta))   draws$zeta  <- zeta
  if (!is.null(shape))  draws$shape <- shape
  if (!is.null(lambda)) draws$lambda <- lambda
  if (!all(weights == 1)) data$weights <- weights
  if (!is.null(offset)) data$offset <- offset
  nlist(data, draws, S = NROW(beta), N = nrow(data))
}

# returns log-likelihood function for loo() and waic()
.llfun <- function(object) {
  f <- object$family
  if (is(f, "family")) get(paste0(".ll_", f$family, "_i"))
  else if (is.character(f)) .ll_polr_i
  else stop("'family' must be a family or a character string")
}

.ll_polr_i <- function(i, data, draws) {
  rmv <- c("y", "weights","offset")
  xdat <- data[, -which(colnames(data) %in% rmv)]
  
  eta <- as.vector(linear_predictor(draws$beta, xdat, data$offset))
  f <- draws$f
  J <- draws$max_y
  y_i <- data$y
  
  if (f == "logistic")    linkinv <- make.link("logit")$linkinv
  else if (f == "loglog") linkinv <- pgumbel
  else                    linkinv <- make.link(f)$linkinv
  
  val <- 
    if (y_i == 1) log(linkinv(draws$zeta[,1] - eta))
    else if (y_i == J) log1p(-linkinv(draws$zeta[,J-1] - eta))
    else log(linkinv(draws$zeta[,y_i] - eta) - 
             linkinv(draws$zeta[,y_i - 1L] - eta))
  
  if ("weights" %in% names(data)) val * data$weights
  else val
}

.ll_gaussian_i <- function(i, data, draws) {
  rmv <- c("y","weights","offset")
  xdat <- data[, -which(colnames(data) %in% rmv)]
  mu <- as.vector(draws$f$linkinv(linear_predictor(draws$beta, xdat, data$offset)))
  val <- dnorm(data$y, mean = mu, sd = draws$sigma, log = TRUE)
  if ("weights" %in% names(data)) val * data$weights
  else val
}
.ll_binomial_i <- function(i, data, draws) {
  rmv <- c("y", "trials", "weights","offset")
  xdat <- data[, -which(colnames(data) %in% rmv)]
  p <- as.vector(draws$f$linkinv(linear_predictor(draws$beta, xdat, data$offset)))
  val <- dbinom(data$y, size = data$trials, prob = p, log = TRUE)
  if ("weights" %in% names(data)) val * data$weights
  else val
}
.ll_poisson_i <- function(i, data, draws) {
  rmv <- c("y", "weights","offset")
  xdat <- data[, -which(colnames(data) %in% rmv)]
  lambda <- as.vector(draws$f$linkinv(linear_predictor(draws$beta, xdat, data$offset)))
  val <- dpois(data$y, lambda, log = TRUE)
  if ("weights" %in% names(data)) val * data$weights
  else val
}
.ll_neg_binomial_2_i <- function(i, data, draws) {
  stop("write the .ll_neg_binomial_2_i function")
}
.ll_Gamma_i <- function(i, data, draws) {
  rmv <- c("y", "weights","offset")
  xdat <- data[, -which(colnames(data) %in% rmv)]
  mu <- as.vector(draws$f$linkinv(linear_predictor(draws$beta, xdat, data$offset)))
  val <- dgamma(data$y, shape = draws$shape, rate = draws$shape / mu, log = TRUE)
  if ("weights" %in% names(data)) val * data$weights
  else val
}
.ll_inverse.gaussian_i <- function(i, data, draws) {
  rmv <- c("y", "weights","offset")
  xdat <- data[, -which(colnames(data) %in% rmv)]
  mu <- as.vector(draws$f$linkinv(linear_predictor(draws$beta, xdat, data$offset)))
  
  val <- 0.5 * log(draws$lambda / (2 * pi)) - 
         1.5 * log(data$y) -
         0.5 * lambda * (data$y - mu)^2 / 
                        (data$y * mu^2)
  if ("weights" %in% names(data)) val * data$weights
  else val
}
