#' Draw from posterior predictive distribution
#' 
#' @export
#' @param object A fitted model object returned by one of the modeling 
#'   functions in this package. This will typically be a list with class 
#'   'stanreg' as well as at least one of 'lm', 'glm', 'polr', or 'lmerMod'.
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used.
#' @param draws The number of draws to return. The default and maximum number of
#'   draws is the size of the posterior sample.
#' @param fun Optional function to apply to the results. See Examples. 
#' 
#' @return A matrix of draws from the posterior predictive distribution.
#' 
#' @examples 
#' fit <- stan_glm(mpg ~ wt, data = mtcars)
#' ppd <- posterior_predict(fit) 
#' mean_wt <- mean(mtcars$wt)
#' hist(posterior_predict(fit, newdata = data.frame(wt = mean_wt)))
#'  
#' fit <- stan_glm(I(log(mpg)) ~ wt, data = mtcars)
#' ppd <- posterior_predict(fit, fun = exp)
#' 
posterior_predict <- function(object, newdata = NULL, draws = NULL, fun) {
  if (object$algorithm == "optimizing")
    stop("posterior_predict only available for MCMC")
  family <- object$family
  famname <- family$family
  ppfun <- paste0(".pp_", famname)
  dat <- .pp_data(object, newdata)
  stanmat <- as.matrix(object$stanfit)
  S <- nrow(stanmat)
  if (is.null(draws)) 
    draws <- S
  else {
    if (draws > S)
      stop(paste("draws =", draws, "but only", S, "draws found."), call. = FALSE)
  } 
  beta <- stanmat[, 1:ncol(dat$x)]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  if (draws < S)
    eta <- eta[sample(S, draws), , drop = FALSE]
  if (is(object, "polr")) {
    f <- object$method
    if (f == "logistic")    linkinv <- make.link("logit")$linkinv
    else if (f == "loglog") linkinv <- pgumbel
    else                    linkinv <- make.link(f)$linkinv
    zeta <- stanmat[,grep("|", colnames(stanmat), value = TRUE, fixed = TRUE)]
    .pp_polr(eta, zeta, linkinv)
  }
  else {
    ppargs <- list(mu = family$linkinv(eta))
    if (famname == "gaussian")
      ppargs$sigma <- stanmat[, "sigma"]
    if (famname == "binomial") {
      y <- if (!is.null(object$y)) 
        object$y else model.response(model.frame(object))
      ppargs$trials <- if (NCOL(y) == 2L) rowSums(y) else rep(1, NROW(y))
    }
    ytilde <- do.call(ppfun, ppargs)
    if (missing(fun)) ytilde
    else do.call(fun, list(ytilde))
  }
}

.pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s,], sigma[s])
  }))
}
.pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s,])
  }))
}
.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s,])
  }))
}
.pp_polr <- function(eta, zeta, linkinv) {
  n <- ncol(eta)
  q <- ncol(zeta)
  t(sapply(1:nrow(eta), FUN = function(s) {
    cumpr <- matrix(linkinv(matrix(zeta[s,], n, q, byrow = TRUE) - eta[s,]), , q)
    fitted <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
    apply(fitted, 1, function(p) which(rmultinom(1, 1, p) == 1))
  }))
}