#' @method print stanreg
#' @export
print.stanreg <- function(x, digits = 3, ...) {
  if (x$algorithm != "optimizing") {
    nms <- setdiff(rownames(x$stan_summary), "log-posterior")
    mat <- as.matrix(x$stanfit)[,nms,drop=FALSE]
    print(cbind(Median = apply(mat, 2, median), MAD_SD = apply(mat, 2, mad)), 
          digits = digits, ...)
  }
  else {
    nms <- names(x$coefficients)
    famname <- x$family$family
    if (is(x, "polr")) 
      nms <- c(nms, grep("|", rownames(x$stan_summary), 
                         fixed = TRUE, value = TRUE))
    else if (is.gaussian(famname))
      nms <- c(nms, "sigma")
    else if (is.gamma(famname))
      nms <- c(nms, "shape")
    else if (is.ig(famname))
      nms <- c(nms, "lambda")
    else if (is.nb(famname))
      nms <- c(nms, "overdispersion")
    nms <- c(nms, grep("^mean_PPD", rownames(x$stan_summary), value = TRUE))
    print(x$stan_summary[nms,1:2], digits = digits, ...)
  }

  if (is(x, "aov")) {
    labels <- attributes(x$terms)$term.labels
    patterns <- gsub(":", ".*:", labels)
    groups <- sapply(patterns, simplify = FALSE, FUN = grep, 
                     x = dimnames(extract(x$stanfit, pars = "beta", permuted = FALSE))$parameters)
    names(groups) <- gsub(".*", "", names(groups), fixed = TRUE)
    groups <- groups[sapply(groups, length) > 0]
    effects_dim <- dim(x$effects)
    effects <- x$effects^2
    effects <- sapply(groups, FUN = function(i) apply(effects[,,i,drop = FALSE], 1:2, mean))
    dim(effects) <- c(effects_dim[-3], ncol(effects))
    dim(effects) <- c(nrow(effects) * ncol(effects), dim(effects)[3])
    colnames(effects) <- paste("Mean Sq", names(groups))
    cat("ANOVA-like table\n")
    print(cbind(Median = apply(effects, 2, median),
                MAD_SD = apply(effects, 2, mad)), digits = digits, ...)
  }
  return(invisible(NULL))
}
