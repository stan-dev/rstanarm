#' Draw from posterior predictive distribution
#' 
#' The posterior predictive distribution is the distribution of the outcome 
#' implied by the model after using the observed data to update our beliefs 
#' about the unknown parameters. Simulating data from the posterior predictive 
#' distribution using the observed predictors is useful for checking the fit of
#' the model. Drawing from the posterior predictive distribution at interesting
#' values of the predictors also lets us visualize how a manipulation of a
#' predictor affects (a function of) the outcome(s).
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used.
#' @param draws The number of draws to return. The default and maximum number of
#'   draws is the size of the posterior sample.
#' @param fun An optional function to apply to the results. This can be the name
#'   of a function as a character string (e.g. \code{fun = "exp"}) or a function
#'   object (e.g. \code{fun = exp}, or \code{fun = function(x) exp(x)}, etc.).
#'   See Examples.
#' @param allow.new.levels A logical scalar (defaulting to \code{FALSE})
#'   indicating whether new levels in grouping variables are allowed in
#'   \code{newdata}. Only relevant for \code{\link[=stan_glmer]{GLMS with
#'   group-specific terms}}.
#' @param ... Currently ignored.
#' 
#' @return A \code{draws} by \code{nrow(newdata)} matrix of simulations
#'   from the posterior predictive distribution. Each row of the matrix is a
#'   vector of predictions generated using a single draw of the model parameters
#'   from the posterior distribution.
#' 
#' @seealso \code{\link{ppcheck}} for graphical posterior predictive checks.
#'   Examples of posterior predictive checking can also be found in the
#'   \pkg{rstanarm} vignettes and demos.
#'   
#' @examples
#' yrep <- posterior_predict(example_model)
#' table(yrep)
#' 
#' \dontrun{
#' nd <- lme4::cbpp
#' nd$size <- max(nd$size) + 1L
#' ppd <- posterior_predict(example_model, newdata = nd)
#' 
#' # Use fun argument to transform predictions
#' fit <- stan_glm(I(log(mpg)) ~ wt, data = mtcars)
#' ppd <- posterior_predict(fit, fun = exp)
#' }
#' 
posterior_predict <- function(object, newdata = NULL, draws = NULL, 
                              allow.new.levels = FALSE, fun, ...) {
  if (!is.stanreg(object))
    stop(deparse(substitute(object)), " is not a stanreg object")
  if (used.optimizing(object)) 
    STOP_not_optimizing("posterior_predict")
  family <- object$family
  if (!is(object, "polr")) {
    famname <- family$family
    ppfun <- paste0(".pp_", famname) 
  }

  S <- .posterior_sample_size(object)
  if (is.null(draws)) draws <- S
  if (draws > S) stop(paste("draws =", draws, "but only", S, "draws found."), 
                      call. = FALSE)
  if (!is.null(newdata)) newdata <- as.data.frame(newdata)
  dat <- pp_data(object, newdata, allow.new.levels, ...)
  x <- dat$x
  if (is.null(NEW_ids <- attr(x, "NEW_ids"))) {
    stanmat <- as.matrix(object)
    beta <- stanmat[, 1:ncol(x), drop = FALSE]
  }
  else { # newdata has new levels
    stanmat <- as.matrix(object$stanfit)
    mark <- grepl("_NEW_", colnames(stanmat), fixed = TRUE)
    if (any(mark)) {
      NEW_draws <- stanmat[, mark, drop = FALSE]
      stanmat <- stanmat[, !mark, drop = FALSE] 
    }
    NEW_cols <- unlist(NEW_ids, use.names = FALSE, recursive = TRUE)
    xNEW <- x[, NEW_cols, drop = FALSE]
    x <- x[, -NEW_cols, drop = FALSE]
    beta <- stanmat[, 1:ncol(x), drop = FALSE]
    x <- cbind(x, xNEW)
    sel <- vector("list", length(NEW_ids))
    for (j in seq_along(NEW_ids)) {
      sel[[j]] <- list()
      idj <- NEW_ids[[j]]
      for (k in seq_along(idj)) {
        patt <- paste0("b[", names(idj)[k],":_NEW_]")
        sel[[j]][[k]] <- rep(grep(patt, colnames(NEW_draws), fixed = TRUE), 
                             length = length(idj[[k]])) # replicate to select same column of NEW_draws multiple times (if necessary)
      }
    }
    beta <- cbind(beta, NEW_draws[, unlist(sel)])
  }
  eta <- linear_predictor(beta, x, dat$offset)
  inverse_link <- linkinv(object)
  if (draws < S)
    eta <- eta[sample(S, draws),, drop = FALSE]
  if (is(object, "polr")) {
    zeta <- stanmat[, grep("|", colnames(stanmat), value = TRUE, fixed = TRUE)]
    ytilde <- .pp_polr(eta, zeta, inverse_link)
  }
  else {
    ppargs <- list(mu = inverse_link(eta))
    if (is.gaussian(famname))
      ppargs$sigma <- stanmat[, "sigma"]
    else if (is.binomial(famname)) {
      y <- get_y(object)
      if (NCOL(y) == 2L) ppargs$trials <- rowSums(y)
      else if (!all(y %in% c(0, 1))) ppargs$trials <- object$weights
      else ppargs$trials <- rep(1, NROW(y))
    }
    else if (is.gamma(famname))
      ppargs$scale <- stanmat[,"scale"]
    else if (is.ig(famname))
      ppargs$lambda <- stanmat[,"lambda"]
    else if (is.nb(famname))
      ppargs$size <- stanmat[,"overdispersion"]
    
    ytilde <- do.call(ppfun, ppargs)
  }
  
  if (!missing(newdata) && nrow(newdata) == 1L) ytilde <- t(ytilde)
  if (!missing(fun)) return(do.call(match.fun(fun), list(ytilde)))
  else return(ytilde)
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
.pp_neg_binomial_2 <- function(mu, size) {
  t(sapply(1:nrow(mu), function(s) {
    rnbinom(ncol(mu), size = size[s], mu = mu[s,])
  }))
}
.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s,])
  }))
}
.pp_Gamma <- function(mu, shape) {
  t(sapply(1:nrow(mu), function(s) {
    rgamma(ncol(mu), shape = shape[s], rate = shape[s] / mu[s,])
  }))
}
.pp_inverse.gaussian <- function(mu, lambda) {
  t(sapply(1:nrow(mu), function(s) {
    .rinvGauss(ncol(mu), mu = mu[s,], lambda = lambda[s])
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

.rinvGauss <- function(n, mu, lambda) {
  mu2 <- mu^2
  y <- rnorm(n)^2
  z <- runif(n)
  x <- mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * y^2) ) / (2 * lambda)
  ifelse (z <= (mu / (mu + x)), x, mu2 / x)
}
