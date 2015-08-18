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
  loo.function(.llfun(x), args = .llargs(x), ...)
}

#' @rdname loo.stanreg
#' @export
#' @importFrom loo waic waic.function
#' @note The \code{...} is ignored for \code{waic}.
#' 
waic.stanreg <- function(x, ...) {
  waic.function(.llfun(x), args = .llargs(x))
}

# returns args argument for loo() and waic()
.llargs <- function(x) {
  args <- list()
  args$family <- x$family
  stanmat <- as.matrix(x$stanfit)
  nms <- names(x)
  args$x <- if ("x" %in% nms) x$x else model.matrix(x)
  args$y <- x$y
  args$offset <- x$offset
  args$weights <- x$weights
  if (is(args$family, "family") && args$family$family == "gaussian") 
    args$sigma <- stanmat[, "sigma"]
  if (is.character(args$family)) {
    args$beta <- stanmat[,colnames(args$x),drop = FALSE]
    args$zeta <- stanmat[,grep("|", colnames(stanmat), fixed = TRUE, value = TRUE),drop=FALSE]
  } else {
    args$beta <- stanmat[, 1:ncol(args$x)]
  }
  do.call(".llargs_prep", args)
}
.llargs_prep <- function(family, x, y, weights, beta,
                         # remaining arguments might be applicable
                         sigma = NULL, zeta = NULL, offset = NULL) {
  eta <- linear_predictor(beta, x, offset)
  f <- family
  if (is(f, "family")) {
    draws <- list(mu = f$linkinv(eta))
    if (f$family != "binomial") data <- nlist(y)
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
      data <- nlist(y, trials)
    } 
  }
  else if (is.character(f)) {
    y <- as.integer(y)
    data <- nlist(y, f, max_y = max(y))
    draws <- nlist(eta)
  }
  else stop("'family' must be a family or a character string")
  
  if (!is.null(sigma)) draws$sigma <- sigma
  if (!is.null(zeta))  draws$zeta  <- zeta
  if (!all(weights == 1)) data$weights <- weights
  nlist(data, draws, S = NROW(beta), N = NROW(y))
}

# returns log-likelihoo function for loo() and waic()
.llfun <- function(object) {
  f <- object$family
  if (is(f, "family")) get(paste0(".ll_", f$family, "_i"))
  else if (is.character(f)) .ll_polr_i
  else stop("'family' must be a family or a character string")
}

.ll_polr_i <- function(i, data, draws) {
  f <- data$f
  J <- data$max_y
  y_i <- data$y[i]
  
  if (f == "logistic")    linkinv <- make.link("logit")$linkinv
  else if (f == "loglog") linkinv <- pgumbel
  else                    linkinv <- make.link(f)$linkinv
  
  val <- 
    if (y_i == 1) log(linkinv(draws$zeta[,1] - draws$eta[,i]))
  else if (y_i == J) log1p(-linkinv(draws$zeta[,J-1] - draws$eta[,i]))
  else log(linkinv(draws$zeta[,y_i] - draws$eta[,i]) - 
             linkinv(draws$zeta[,y_i - 1L] - draws$eta[,i]))
  
  if ("weights" %in% names(data)) val * data$weights[i]
  else val
}

.ll_gaussian_i <- function(i, data, draws) {
  val <- dnorm(data$y[i], mean = draws$mu[,i], sd = draws$sigma, log = TRUE)
  if ("weights" %in% names(data)) val * data$weights[i]
  else val
}
.ll_binomial_i <- function(i, data, draws) {
  val <- dbinom(data$y[i], size = data$trials[i], prob = draws$mu[,i], log = TRUE)
  if ("weights" %in% names(data)) val * data$weights[i]
  else val
}
.ll_poisson_i <- function(i, data, draws) {
  val <- dpois(data$y[i], lambda = draws$mu[,i], log = TRUE)
  if ("weights" %in% names(data)) val * data$weights[i]
  else val
}
