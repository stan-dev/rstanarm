#' @export
stan_nlmer <- function (formula, data = NULL, subset, weights, na.action, offset, 
                        contrasts = NULL, ..., 
                        prior = normal(), prior_ops = prior_options(),
                        prior_covariance = decov(), prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  f <- as.character(formula[-3])
  SSfunctions <- grep("^SS[[:lower:]]+", ls("package:stats"), value = TRUE)
  SSfun <- sapply(SSfunctions, grepl, x = f[2])
  if (any(SSfun)) SSfun <- which(SSfun)
  else stop("'stan_nlmer' requires a named self-starting nonlinear function")
  mc <- match.call(expand.dots = FALSE)
  mc$start <- getInitial(as.formula(f), data)
  nlf <- nlformula(mc)
  y <- nlf$respMod$y
  X <- nlf$X
  nlf$reTrms$SSfun <- SSfun
  if (SSfun == 5) {
    nlf$reTrms$Dose <- nlf$frame[[2]]
    nlf$reTrms$input <- nlf$frame[[3]]
  }
  else nlf$reTrms$input <- nlf$frame[[2]]
  
  nlf$reTrms$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, family = gaussian(),
                          weights = double(), offset = double(), # FIXME
                          prior = prior, prior_intercept = NULL,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = nlf$reTrms, QR = QR, sparse = sparse, ...)
  return(stanfit)
}
