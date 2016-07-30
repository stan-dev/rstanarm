args(betareg::betareg)

stan_betareg <- function (formula, data, subset, na.action, weights, offset,
                          link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                          link.phi = NULL, model = TRUE, y = TRUE, x = FALSE, ...,
                          prior = normal(), prior_intercept = normal(),
                          prior_ops = prior_options(), prior_PD = FALSE, 
                          algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                          adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  
  mc <- match.call(expand.dots = FALSE)
  mc$model <- mc$y <- mc$x <- TRUE
  
  # NULLify any Stan specific arguments in mc now
  mc$prior <- mc$prior_intercept <- mc$prior_ops <- mc$prior_PD <- mc$algorithm <-
    mc$adapt_delta <- mc$QR <- mc$sparse <- NULL
  
  mc$drop.unused.levels <- TRUE
  mc[[1L]] <- quote(betareg::betareg)
  
  if (!requireNamespace("betareg")) stop("the betareg package is needed by 'stan_betareg'")
  mc$control <- betareg::betareg.control(maxit = 0, fsmaxit = 0)
  br <- suppressWarnings(eval(mc, parent.frame()))
  mf <- check_constant_vars(br$model)
  mt <- br$terms
  Y <- array1D_check(model.response(mf, type = "any"))
  X <- model.matrix(br)
  Z <- model.matrix(br, model = "precision")
  #if(ncol(Z) == 1 && all(Z == 1)) Z <- NULL
  
  weights <- validate_weights(as.vector(model.weights(mf)))
  offset <- validate_offset(as.vector(model.offset(mf)), y = Y)
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  
  # pass existence of declaration of linear predictor of the dispertion parameter
  Z_true <- all.names(formula)[3]
  if (Z_true=="|") {
    Z_true <- 1
  }
  else {
    Z_true <- 0
  }

  
  # pass the prior information to stan_betareg.fit()
  stanfit <- stan_betareg.fit(x = X, y = Y, z = Z, weights = NULL, offset = NULL, 
                              link = link, link.phi = link.phi, ..., prior = prior, 
                              prior_intercept = prior_intercept, prior_ops = prior_ops,
                              prior_PD = prior_PD, algorithm = algorithm, 
                              adapt_delta = adapt_delta, QR = QR, sparse = FALSE, Z_true = Z_true)
  algorithm <- match.arg(algorithm)
  link <- match.arg(link)
  fit <- nlist(stanfit, family = beta_fam(link), formula, offset = NULL, 
               weights = NULL, x = X, y = Y, 
               data, prior.info = get_prior_info(call, formals()), 
               call = match.call(), terms = mt, model = mf, 
               algorithm, na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"))
  out <- stanreg(fit)
  class(out) <- c("stanreg", "betareg")
  # out$xlevels <- .getXlevels(mt, mf)
  if (!x) 
    out$x <- NULL
  if (!y) 
    out$y <- NULL
  if (!model) 
    out$model <- NULL
  
  return(out) 
}

beta_fam <- function(link = "logit") {
  stopifnot(is.character(link))
  if (link == "loglog") {
    out <- binomial("cloglog")
    out$linkinv <- function(eta) 1 - pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
    out$linkfun <- function(mu) log(-log(mu))
  }
  else out <- binomial(link)
  out$family <- "beta"
  out$variance <- function(mu, phi) mu * (1 - mu) / (phi + 1)
  out$dev.resids <- function(y, mu, wt)
    stop("'dev.resids' function should not be called")
  out$aic <- function(y, n, mu, wt, dev)
    stop("'aic' function should not have been called")
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}
