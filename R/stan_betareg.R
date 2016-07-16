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
  mc[[1L]] <- quote(betareg::betareg)
  if (!requireNamespace("betareg")) stop("the betareg package is needed by 'stan_betareg'")
  mc$control <- betareg::betareg.control(maxit = 0, fsmaxit = 0)
  br <- suppressWarnings(eval(mc, parent.frame()))
  mf <- check_constant_vars(br$model)
  mt <- attr(mf, "terms")
  Y <- array1D_check(model.response(mf, type = "any"))
  if (is.empty.model(mt))
    stop("No intercept or predictors specified.", call. = FALSE)
  X <- model.matrix(mt, mf, contrasts)
  
  # pass the prior information to stan_betareg.fit()
  stanfit <- stan_betareg.fit(x = X, y = Y, z = NULL, weights = NULL, offset = NULL, 
                              link = link, link.phi = link.phi, ..., prior = prior, 
                              prior_intercept = prior_intercept, prior_ops = prior_ops,
                              prior_PD = prior_PD, algorithm = algorithm, 
                              adapt_delta = adapt_delta, QR = QR, sparse = FALSE)
  algorithm <- match.arg(algorithm)
  link <- match.arg(link)
  fam <- binomial(link = link)
  fam$family <- "beta"
  fit <- nlist(stanfit, family = fam, formula, offset = NULL, 
               weights = NULL, x = X, y = Y, 
               data, prior.info = get_prior_info(call, formals()), 
               call = match.call(), terms = mt, model = mf, 
               algorithm, na.action = attr(mf, "na.action"), 
               contrasts = attr(X, "contrasts"))
  out <- stanreg(fit)
  out$xlevels <- .getXlevels(mt, mf)
  if (!x) 
    out$x <- NULL
  if (!y) 
    out$y <- NULL
  if (!model) 
    out$model <- NULL
  
  return(out) 
}
