.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

#' @method print stanreg
#' @export
print.stanreg <- function(x, digits = 1, ...) {
  print(x$call)
  cat("\nEstimates:\n")
  
  mer <- is(x, "lmerMod")
  ord <- is(x, "polr")
  if (x$algorithm != "optimizing") {
    mat <- as.matrix(x$stanfit)
    nms <- setdiff(rownames(x$stan_summary), "log-posterior")
    if (mer) nms <- setdiff(nms, grep("^b\\[", nms, value = TRUE))
    if (ord) {
      cut_nms <- grep("|", nms, fixed = TRUE, value = TRUE)
      nms <- setdiff(nms, cut_nms)
      cut_mat <- mat[, cut_nms, drop = FALSE]
      cut_estimates <- .median_and_madsd(cut_mat)
    }
    ppd_nms <- grep("^mean_PPD", nms, value = TRUE)
    nms <- setdiff(nms, ppd_nms)
    coef_mat <- mat[, nms, drop = FALSE]
    ppd_mat <- mat[, ppd_nms, drop = FALSE]
    estimates <- .median_and_madsd(coef_mat)
    ppd_estimates <- .median_and_madsd(ppd_mat)
    
    .printfr(estimates, digits, ...)
    if (ord) {
      cat("\nCutpoints:\n")
      .printfr(cut_estimates, digits, ...)
    }
    if (mer) {
      cat("\nError terms:\n")
      print(VarCorr(x), digits = digits + 1, ...)
      cat("Num. levels:", 
          paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
    }
    cat("\nSample avg. posterior predictive \ndistribution of y (X = xbar):\n")
    .printfr(ppd_estimates, digits, ...)
  }
  else {
    nms <- names(x$coefficients)
    if (ord) 
      nms <- c(nms, grep("|", rownames(x$stan_summary), 
                         fixed = TRUE, value = TRUE))
    else {
      famname <- x$family$family
      if (is.gaussian(famname))
        nms <- c(nms, "sigma")
      else if (is.gamma(famname))
        nms <- c(nms, "shape")
      else if (is.ig(famname))
        nms <- c(nms, "lambda")
      else if (is.nb(famname))
        nms <- c(nms, "overdispersion")
    }
    nms <- c(nms, grep("^mean_PPD", rownames(x$stan_summary), value = TRUE))
    estimates <- x$stan_summary[nms,1:2]
    .printfr(estimates, digits, ...)
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
    cat("\nANOVA-like table:\n")
    anova_table <- .median_and_madsd(effects)
    .printfr(anova_table, digits, ...)
  }
  return(invisible(NULL))
}
