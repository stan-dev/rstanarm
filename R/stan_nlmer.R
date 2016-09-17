#' @export
stan_nlmer <- function (formula, data = NULL, subset, weights, na.action, offset, 
                        contrasts = NULL, ..., 
                        prior = normal(), prior_ops = prior_options(),
                        prior_covariance = decov(), prior_PD = FALSE, 
                        algorithm = c("sampling", "meanfield", "fullrank"), 
                        adapt_delta = NULL, sparse = FALSE) {
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
  inputs <- as.character(nlf$respMod$nlmod[2])
  inputs <- sub("(", ",", inputs, fixed = TRUE)
  inputs <- sub(")", "", inputs, fixed = TRUE)
  inputs <- scan(text = inputs, what = character(), sep = ",", strip.white = TRUE, quiet = TRUE)
  if (SSfun == 5) {
    nlf$reTrms$Dose <- nlf$frame[[inputs[2]]]
    nlf$reTrms$input <- nlf$frame[[inputs[3]]]
  }
  else nlf$reTrms$input <- nlf$frame[[inputs[2]]]
  
  nlf$reTrms$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  g <- gaussian(link = "identity")
  weights <- nlf$respMod$weights
  offset  <- nlf$respMod$offset
  stanfit <- stan_glm.fit(x = X, y = y, family = g,
                          weights = weights, offset = offset,
                          prior = prior, prior_intercept = NULL,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = nlf$reTrms, QR = FALSE, sparse = sparse, ...)
  
  Z <- pad_reTrms(Z = t(nlf$reTrms$Zt), cnms = nlf$reTrms$cnms, 
                  flist = nlf$reTrms$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  g$link <- paste("inv", SSfunctions[SSfun], sep = "_")
  g$linkinv <- function(eta) {
    SSargs <- as.data.frame(matrix(eta, nrow = length(y), ncol = ncol(X), 
                                   dimnames = list(NULL, colnames(X))))
    if (SSfun == 5) SSargs <- cbind(Dose = nlf$frame[[inputs[2]]], 
                                    input = nlf$frame[[inputs[3]]], SSargs)
    else if (SSfun == 7 || SSfun == 10) SSargs <- cbind(x = nlf$frame[[inputs[2]]], SSargs)
    else SSargs <- cbind(input = nlf$frame[[inputs[2]]], SSargs)
    do.call(SSfunctions[SSfun], args = SSargs)
  }
  g$linkfun  <- function(mu) stop("'linkfun' should not have been called")
  g$variance <- function(mu) stop("'variance' should not have been called")
  g$mu.eta   <- function(mu) stop("'mu.eta' should not have been called")
  fit <- nlist(stanfit, family = g, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call = match.call(), terms = NULL, model = NULL, 
               prior.info = get_prior_info(call, formals()),
               na.action = na.omit, contrasts, algorithm, glmod = nlf)
  out <- stanreg(fit)
  class(out) <- c(class(out), "nlmerMod", "lmerMod")
  return(out)
}
