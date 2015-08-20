#' @method summary stanreg
#' @export
summary.stanreg <- function(object, ..., digits = 2) {
  # use RStan's summary just as placeholder. we should replace this with our own
  # summary 
  if (object$stanfit@mode == 0) {
    out <- round(summary(object$stanfit, ...)$summary, digits = digits)
    if ("n_eff" %in% colnames(out))
      out[, "n_eff"] <- round(out[, "n_eff"])
    if ("se_mean" %in% colnames(out)) # So people don't confuse se_mean and sd
      colnames(out)[which(colnames(out) == "se_mean")] <- "mcse" 
    out
  }
  else {
    mark <- names(object$coefficients)
    if (object$family$family == "gaussian") mark <- c(mark, "sigma")
    else if (object$family$family == "Negative Binomial") mark <- c(mark, "overdispersion")
    object$stan_summary[mark,,drop=FALSE]
  }
}
