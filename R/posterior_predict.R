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
#'   transformations were specified inside the model formula.
#' @param draws The number of draws to return. The default and maximum number of
#'   draws is the size of the posterior sample.
#' @param fun An optional function to apply to the results. \code{fun} is found 
#'   by a call to \code{\link{match.fun}} and so can be specified as a function
#'   object, a string naming a function, etc.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently, \code{...} can contain two arguments that are only 
#' relevant for \code{\link[=stan_glmer]{GLMS with
#' group-specific terms}}:  
#' \describe{
#' \item{\code{allow.new.levels}}{A logical scalar (defaulting to \code{FALSE}) 
#' indicating whether new levels in grouping variables are allowed in 
#' \code{newdata}.}
#' \item{\code{re.form}}{A formula for "random effects" to condition on. If
#' \code{NULL} (the default), all are included. If \code{NA} or \code{~0}, all
#' are omitted.}
#'}
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
                              fun, seed, ...) {
  if (!missing(seed)) 
    set.seed(seed)
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
  if (draws > S) {
    stop(paste0("'draws' = ", draws, 
                " but posterior sample size is only ", S, "."))
  }
  if (!is.null(newdata)) {
    newdata <- as.data.frame(newdata)
    if (any(is.na(newdata))) 
      stop("Currently NAs are not allowed in 'newdata'.")
  }
  dat <- pp_data(object, newdata, ...)
  x <- dat$x
  if (is.null(attr(x, "NEW_ids"))) { # no new levels in grouping variables
    stanmat <- as.matrix(object)
    beta <- stanmat[, 1:ncol(x), drop = FALSE]
  } else { # newdata has new levels
    stanmat <- as.matrix(object$stanfit)
    tmp <- pp_new_levels(stanmat, x)
    x <- tmp$x
    beta <- tmp$beta
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
      ppargs$shape <- stanmat[,"shape"]
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
.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s,])
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

# Create x and beta when new grouping levels are present 
# 
# @param stanmat Matrix of posterior draws, from calling
#   as.matrix(object$stanfit) instead of as.matrix(object) since we need to
#   include the extra "_NEW_" parameters used for predicting new levels.
# @param x The x matrix returned by \code{.pp_data_mer} (which also includes the
#   Z matrix). \code{x} should have a "NEW_ids" attribute.
# @return A named list containing \code{x} and \code{beta}.
pp_new_levels <- function(stanmat, x) {
  NEW_ids <- attr(x, "NEW_ids")
  mark <- grepl("_NEW_", colnames(stanmat), fixed = TRUE)
  if (is.null(NEW_ids)) 
    stop("Missing 'NEW_ids' attribute necessary for prediction.")
  if (!any(mark)) 
    stop("Draws used to predict for new levels were not found.")
  
  NEW_draws <- stanmat[, mark, drop = FALSE]
  stanmat <- stanmat[, !mark, drop = FALSE] 
  NEW_cols <- unlist(NEW_ids, use.names = FALSE, recursive = TRUE)
  xNEW <- x[, NEW_cols, drop = FALSE]
  x <- x[, -NEW_cols, drop = FALSE]
  beta <- stanmat[, 1:ncol(x), drop = FALSE]
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
  betaNEW <- NEW_draws[, unlist(sel), drop = FALSE]
  
  nc <- ncol(x) + ncol(xNEW)
  mark <- c(OLD_cols = setdiff(1:nc, NEW_cols), NEW_cols)
  xout <- matrix(NA, nrow = nrow(x), ncol = nc)
  bout <- matrix(NA, nrow = nrow(beta), ncol = nc)
  xout[, mark] <- cbind(x, xNEW)
  bout[, mark] <- cbind(beta, betaNEW)
  colnames(xout)[mark] <- c(colnames(x), colnames(xNEW))
  colnames(bout)[mark] <- c(colnames(beta), colnames(betaNEW))
  return(list(x = xout, beta = bout))
}

# this is mostly copied from lme4:::predict.merMod and lme4:::mkNewReTrms
new_thing <- function(object, newdata = NULL, draws = NULL, 
                      fun, seed, re.form, ...) {
  X <- object$x[,1:length(fixef(object)), drop = FALSE]
  R <- formula(object, fixed.only = TRUE)
  R <- R[[length(R)]]
  RHS <- formula(substitute(~R, list(R = R)))
  Terms <- terms(object, fixed.only = TRUE)
  mf <- model.frame(object, fixed.only = TRUE)
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) NULL else lapply(mf[isFac], levels)
  mfnew <- model.frame(delete.response(Terms), newdata, xlev = orig_levs)
  X <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(X, "contrasts"))
  offset <- 0
  tt <- terms(object)
  if (!is.null(off.num <- attr(tt, "offset"))) {
    for (i in off.num)
      offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
  }
  stanmat <- as.matrix(object$stanfit)
  beta <- stanmat[,1:ncol(X), drop = FALSE]
  eta <- linear_predictor(beta, X, offset)
  if (!missing(re.form)) {
    if (is.null(newdata)) rfd <- mfnew <- model.frame(object)
    else {
      mfnew <- model.frame(delete.response(terms(object, fixed.only = TRUE)), 
                           newdata)
      tt <- delete.response(terms(object, random.only = TRUE))
      rfd <- model.frame(tt, newdata, na.action = na.pass)
    }
    if (inherits(re.form, "formula")) {
      ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
      ns.re <- names(re <- ranef(object))
      nRnms <- names(Rcnms <- ReTrms$cnms)
      if (!all(nRnms %in% ns.re)) 
        stop("grouping factors specified in re.form that were not present in original model")
      new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
      Zt <- ReTrms$Zt
      p <- sapply(ReTrms$cnms, FUN = length)
      l <- sapply(attr(ReTrms$flist, "assign"), function(i) 
        nlevels(ReTrms$flist[[i]]))
      t <- length(p)
      group_nms <- names(ReTrms$cnms)
      Z_names <- unlist(lapply(1:t, FUN = function(i) {
        paste0(ReTrms$cnms[[i]], " ", group_nms[i], ":", levels(ReTrms$flist[[i]]))
      }))
      b <- stanmat[,grepl("^b\\[", colnames(stanmat)), drop = FALSE]
      ord <- sapply(Z_names, FUN = function(x) {
        m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
        len <- length(m)
        if (len == 1) return(m)
        if (len > 1) stop("multiple matches bug")
        x <- sub(":.*$", ":_NEW_", x)
        grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
      })
      b <- b[,ord, drop = FALSE]
      eta <- eta + as.matrix(b %*% Zt)
    }
  }
  inverse_link <- linkinv(object)
  family <- object$family
  if (!is(object, "polr")) {
    famname <- family$family
    ppfun <- paste0(".pp_", famname) 
  }
  S <- .posterior_sample_size(object)
  if (is.null(draws)) draws <- S
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
      ppargs$shape <- stanmat[,"shape"]
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
