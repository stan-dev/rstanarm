stan_nlmer <- function (formula, data = NULL, subset, weights, na.action, offset, 
                        contrasts = NULL, ..., 
                        prior = normal(), prior_ops = prior_options(),
                        prior_covariance = decov(), prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta = NULL, QR = FALSE, sparse = FALSE) {
  f <- as.character(formula[-3])
  SSfunctions <- grep("^SS[[:lower:]]+", ls("package:stats"), value = TRUE)
  SSfun <- sapply(SSfunctions, grepl, x = f)
  if (any(SSfun)) SSfun <- which(SSfun)
  else stop("'stan_nlmer' requires a named self-starting nonlinear function")
  mc <- match.call(expand.dots = FALSE)
  mc$start <- getInitial(as.formula(f), data)
  nlf <- nlformula(mc)
  y <- nlf$frame[[1]]
  X <- nlf$X
  Z <- t(nlf$reTrms$Zt)
  nlf$reTrms$SSfun <- SSfun
  return(nlf$reTrms)
}
