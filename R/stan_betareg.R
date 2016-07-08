args(betareg::betareg)

stan_betareg <- function (formula, data, subset, na.action, weights, offset, 
          link = c("logit", "probit", "cloglog", "cauchit", "log", 
                   "loglog"), link.phi = NULL, model = TRUE, y = TRUE, x = FALSE, ...,
          prior = normal(), prior_intercept = normal(),
          prior_ops = prior_options(), prior_PD = FALSE, 
          algorithm = c("sampling", "optimizing", 
                        "meanfield", "fullrank"),
          adapt_delta = NULL, QR = FALSE, sparse = FALSE)
          )
  {
  mc <- match.call(expand.dots = FALSE)
  mc$model <- mc$y <- mc$x <- TRUE
  # NULLify any Stan specific arguments in mc now
  mc$prior <- mc$prior_intercept <- mc$prior_ops <- mc$prior_PD <- mc$algorithm <-
    mc$adapt_delta <- mc$QR <- mc$sparse <- NULL
  mc[[1L]] <- quote(betareg::betareg)
  if (!requireNamespace("betareg")) stop("the betareg package is needed by 'stan_betareg'")
  mc$betareg.control <- betareg::betareg.control(maxit = 0, fsmaxit = 0)
  mf <- eval(mc, parent.frame())
  
  # pass the prior information to stan_betareg.fit()
  stanfit <- stan_betareg.fit()
  fit <- nlist(stanfit, family = "beta", formula, offset, weights, x = X, y = Y, 
               data, prior.info = get_prior_info(call, formals()), 
               call = call, terms = mt, model = mf, 
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
