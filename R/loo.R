#' Leave-one-out cross-validation (LOO)
#' 
#' Compute approximate leave-one-out cross-validation (LOO) or the Widely 
#' Applicable Information Criterion (WAIC) using the
#' \pkg{\link[=loo-package]{loo}} package.
#' 
#' @export
#' @aliases loo waic
#' @inheritParams loo::loo
#' @param x A fitted model object returned by one of the \pkg{rstanarm} modeling
#'   functions.
#' @return An object of class 'loo'. See \code{\link[loo]{loo}} and
#'   \code{\link[loo]{waic}}.
#' 
#' @examples 
#' seed <- 42024
#' set.seed(seed)
#' fit1 <- stan_glm(mpg ~ wt, data = mtcars, iter = 200, seed = seed)
#' fit2 <- stan_glm(mpg ~ ., data = mtcars, iter = 200, seed = seed)
#' loo1 <- loo(fit1)
#' loo2 <- loo(fit2)
#' loo::compare(loo1, loo2)
#' plot(loo2, label_points = TRUE)
#' 
#' @seealso \code{\link[loo]{loo-package}}
#' @importFrom loo loo loo.function
#' 
loo.stanreg <- function(x, ...) {
  if (x$algorithm != "sampling")
    stop("Only available for MCMC (algorithm = 'sampling').", call. = FALSE)
  loo.function(.llfun(x$family), args = .llargs(x), ...)
}

#' @rdname loo.stanreg
#' @export
#' @importFrom loo waic waic.function
#' @note The \code{...} is ignored for \code{waic}.
#' 
waic.stanreg <- function(x, ...) {
  if (x$algorithm != "sampling")
    stop("Only available for MCMC (algorithm = 'sampling').", call. = FALSE)
  waic.function(.llfun(x$family), args = .llargs(x))
}


# returns log-likelihood function for loo() and waic()
.llfun <- function(f) {
  # f <- object$family
  if (is(f, "family")) get(paste0(".ll_", f$family, "_i"))
  else if (is.character(f)) .ll_polr_i
  else stop("'family' must be a family or a character string")
}

# returns args argument for loo() and waic()
.llargs <- function(object) {
  f <- object$family
  draws <- nlist(f)
  stanmat <- as.matrix(object$stanfit)
  nms <- names(object)
  x <- if ("x" %in% nms) object$x else model.matrix(object)
  y <- if ("y" %in% nms) object$y else model.response(model.frame(object))
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
    draws$beta <- stanmat[, 1:ncol(x)]
    if (f$family == "gaussian") draws$sigma <- stanmat[, "sigma"]
    if (f$family == "Gamma") draws$shape <- stanmat[, "shape"]
    if (f$family == "inverse.gaussian") draws$lambda <- stanmat[, "lambda"]
  }
  else if (is.character(f)) {
    stopifnot(is(object, "polr"))
    y <- as.integer(y)
    data <- data.frame(y, x)
    draws$beta <- stanmat[, colnames(x), drop = FALSE]
    draws$zeta <- stanmat[,grep("|", colnames(stanmat), fixed = TRUE, value = TRUE),drop=FALSE]
    draws$max_y <- max(y)
  }
  else stop("'family' must be a family or a character string")
  
  data$offset <- object$offset
  if (!all(object$weights == 1)) data$weights <- object$weights
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
  stop("write the .ll_neg_binomial_2_i function")
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
  
  if (f == "logistic") linkinv <- make.link("logit")$linkinv
  else if (f == "loglog") linkinv <- pgumbel
  else linkinv <- make.link(f)$linkinv
  
  val <- 
    if (y_i == 1) log(linkinv(draws$zeta[,1] - eta))
    else if (y_i == J) log1p(-linkinv(draws$zeta[,J-1] - eta))
    else log(linkinv(draws$zeta[,y_i] - eta) - 
               linkinv(draws$zeta[,y_i - 1L] - eta))
  
  .weighted(val, data$weights)
}
