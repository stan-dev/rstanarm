# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2017, 2018 Sam Brilleman
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Print method for stanmvreg objects
#
#' @rdname print.stanreg
#' @export
#' @method print stanmvreg
print.stanmvreg <- function(x, digits = 3, ...) {
  
  M     <- get_M(x)
  mvmer <- is.mvmer(x)
  surv  <- is.surv(x)
  jm    <- is.jm(x)
  stubs <- paste0("(", get_stub(x), 1:M, "):")
  
  cat(x$stan_function)
  
  if (mvmer) {
    for (m in 1:M) {
      cat("\n formula", stubs[m], formula_string(formula(x, m = m)))
      cat("\n family ", stubs[m], family_plus_link(x, m = m))
    }    
  }
  
  if (surv) {
    cat("\n formula (Event):", formula_string(formula(x, m = "Event")))
    cat("\n baseline hazard:", get_basehaz_name(x)) 
  }
  
  if (jm) {
    cat("\n assoc:          ", assoc_string(x))
  }
  
  cat("\n------\n")
  
  mat    <- as.matrix(x$stanfit)
  rownms <- rownames(x$stan_summary)
  nms    <- collect_nms(rownms, M, stub = get_stub(x), value = TRUE)
  
  # Estimates table for longitudinal submodel(s)
  if (mvmer) {
    link <- sapply(1:M, function(m) x$family[[m]]$link)
    for (m in 1:M) {
      terms_m <- terms(x)[[m]]
      sel <- attr(terms_m, "response")
      yvar <- rownames(attr(terms_m, "factors"))[sel]
      if (is.jm(x)) {
        cat(paste0("\nLongitudinal submodel", if (M > 1) paste0(" ", m), 
                   ": ", yvar,"\n"))
      } else {
        cat(paste0("\nSubmodel for y", m, ": ", yvar,"\n"))
      }
      coef_mat <- mat[, c(nms$y[[m]], nms$y_extra[[m]]), drop = FALSE]
      
      # Calculate median and MAD
      estimates <- .median_and_madsd(coef_mat)
      
      # Add column with eform
      if (link[m] %in% c("log", "logit")) 
        estimates <- cbind(estimates, 
                           "exp(Median)" = c(exp(estimates[nms$y[[m]], "Median"]), 
                                             rep(NA, length(nms$y_extra[[m]]))))
      
      # Print estimates
      rownames(estimates) <- 
        gsub(paste0("^", get_stub(x), m, "\\|"), "", rownames(estimates))     
      .printfr(estimates, digits, ...)
    }    
  }
  
  # Estimates table for event submodel
  if (surv) {
    cat("\nEvent submodel:\n")   
    coef_mat <- mat[, c(nms$e, nms$a, nms$e_extra), drop = FALSE]
    
    # Calculate median and MAD
    estimates <- .median_and_madsd(coef_mat)
    
    # Add column with eform
    estimates <- cbind(estimates, 
                       "exp(Median)" = c(exp(estimates[c(nms$e, nms$a), "Median"]), 
                                         rep(NA, length(nms$e_extra))))
    
    rownames(estimates) <- gsub("^Event\\|", "", rownames(estimates))  
    rownames(estimates) <- gsub("^Assoc\\|", "", rownames(estimates))   
    .printfr(estimates, digits, ...)
  }
  
  # Estimates table for group-level random effects
  if (mvmer) {
    cat("\nGroup-level error terms:\n") 
    print(VarCorr(x), digits = digits + 1, ...)
    cat("Num. levels:", paste(names(ngrps(x)), unname(ngrps(x)), 
                              collapse = ", "), "\n")  
    
    # Sample average of the PPD
    ppd_mat <- mat[, nms$ppd, drop = FALSE]
    ppd_estimates <- .median_and_madsd(ppd_mat)
    cat("\nSample avg. posterior predictive distribution \nof",
        if (is.jm(x)) "longitudinal outcomes:\n" else "y:\n")
    .printfr(ppd_estimates, digits, ...)
  }
  
  cat("\n------")
  cat("\n* For help interpreting the printed output see ?print.stanreg")
  cat("\n* For info on the priors used see ?prior_summary.stanreg")
  
  invisible(x)
}

#' @rdname summary.stanreg
#' @export
#' @method summary stanmvreg
summary.stanmvreg <- function(object, pars = NULL, regex_pars = NULL, 
                              probs = NULL, ..., digits = 3) {
  
  pars  <- collect_pars(object, pars, regex_pars)
  M     <- get_M(object)
  mvmer <- is.mvmer(object)
  surv  <- is.surv(object)
  jm    <- is.jm(object)  
  
  if (mvmer) {
    # Outcome variable for each longitudinal submodel
    y_vars <- sapply(1:M, function(m, object) {
      terms_m <- terms(object)[[m]]
      sel <- attr(terms_m, "response")
      ret <- rownames(attr(terms_m, "factors"))[sel]
    }, object = object)
    
    # Family and link for each longitudinal submodel
    fam <- lapply(1:M, function(m) family_plus_link(object, m = m))
  }
  
  if (jm) {
    # Association structure
    sel <- grep("^which", rownames(object$assoc), invert = TRUE, value = TRUE)
    assoc <- list_nms(lapply(1:M, function(m) 
      sel[which(object$assoc[sel,m] == TRUE)]), M) 
  }
  
  # Construct summary table  
  args <- list(object = object$stanfit)
  if (!is.null(probs)) 
    args$probs <- probs
  out <- do.call("summary", args)$summary
  nms <- collect_nms(rownames(object$stan_summary), M, 
                     stub = get_stub(object), value = TRUE)
  if (!is.null(pars)) {
    pars2 <- NA     
    if ("alpha" %in% pars) pars2 <- c(pars2, nms$alpha)
    if ("beta"  %in% pars) pars2 <- c(pars2, nms$beta)
    if ("long"  %in% pars) pars2 <- c(pars2, unlist(nms$y), unlist(nms$y_extra))
    if ("event" %in% pars) pars2 <- c(pars2, nms$e, nms$a, nms$e_extra)
    if ("assoc" %in% pars) pars2 <- c(pars2, nms$a)      
    if ("fixef" %in% pars) pars2 <- c(pars2, unlist(nms$y), nms$e, nms$a)
    if ("b"     %in% pars) pars2 <- c(pars2, nms$b)
    pars2 <- c(pars2, setdiff(pars, 
                              c("alpha", "beta", "varying", "b",
                                "long", "event", "assoc", "fixef")))
    pars <- pars2[!is.na(pars2)]
  } else {
    pars <- rownames(object$stan_summary)
    pars <- setdiff(pars, b_names(pars, value = TRUE))
    if (used.variational(object)) 
      pars <- setdiff(pars, "log-posterior")
  }
  out <- out[rownames(out) %in% pars, , drop = FALSE]
  out <- out[!grepl(":_NEW_", rownames(out), fixed = TRUE), , drop = FALSE]
  stats <- colnames(out)
  if ("n_eff" %in% stats)
    out[, "n_eff"] <- round(out[, "n_eff"])
  if ("se_mean" %in% stats) # So people don't confuse se_mean and sd
    colnames(out)[stats %in% "se_mean"] <- "mcse"
  
  # Reorder rows of output table
  nms_tmp <- rownames(out)  
  nms_tmp_y <- lapply(1:M, function(m) 
    grep(paste0("^", get_stub(object), m, "\\|"), nms_tmp, value = TRUE))
  nms_tmp_e <- grep("^Event\\|", nms_tmp, value = TRUE)
  nms_tmp_a <- grep("^Assoc\\|", nms_tmp, value = TRUE)
  nms_tmp_b <- b_names(nms_tmp, value = TRUE)
  nms_tmp_Sigma <- grep("^Sigma", nms_tmp, value = TRUE)
  nms_tmp_lp <- grep("^log-posterior$", nms_tmp, value = TRUE)
  out <- out[c(unlist(nms_tmp_y), nms_tmp_e, nms_tmp_a, nms_tmp_b, 
               nms_tmp_Sigma, nms_tmp_lp), , drop = FALSE]
  
  # Output object
  if (mvmer) {
    out <- structure(out, 
                     y_vars     = y_vars, 
                     family     = fam, 
                     n_markers  = object$n_markers,
                     n_yobs     = object$n_yobs, 
                     n_grps     = object$n_grps)
  }
  if (surv) {
    out <- structure(out, 
                     n_subjects = object$n_subjects, 
                     n_events   = object$n_events,
                     basehaz    = object$survmod$basehaz) 
  }
  if (jm) {
    out <- structure(out, 
                     id_var     = object$id_var, 
                     time_var   = object$time_var, 
                     assoc      = assoc)
  }
  
  structure(out, 
            formula               = object$formula, 
            algorithm             = object$algorithm,
            stan_function         = object$stan_function,
            posterior_sample_size = posterior_sample_size(object),
            runtime               = object$runtime, print.digits = digits,
            class                 = c("summary.stanmvreg", "summary.stanreg"))
}


# Summary method for stanmvreg objects
#
#' @rdname summary.stanreg
#' @export
#' @method print summary.stanmvreg
print.summary.stanmvreg <- function(x, digits = max(1, attr(x, "print.digits")), 
                                    ...) {
  atts  <- attributes(x)
  mvmer <- atts$stan_function %in% c("stan_mvmer", "stan_jm")
  jm    <- atts$stan_function == "stan_jm"
  tab   <- if (jm) "   " else ""
  
  cat("\nModel Info:\n")
  cat("\n function:   ", tab, atts$stan_function)
  if (mvmer) {
    M <- atts$n_markers
    stubs <- paste0("(", if (jm) "Long" else "y", 1:M, "):") 
    for (m in 1:M) {
      cat("\n formula", stubs[m], formula_string(atts$formula[[m]]))
      cat("\n family ", stubs[m], atts$family[[m]])
    }    
  }
  if (jm) {
    assoc_fmt <- uapply(1:M, function(m) paste0(atts$assoc[[m]], " (Long", m, ")"))
    cat("\n formula (Event):", formula_string(atts$formula[["Event"]]))
    cat("\n baseline hazard:", get_basehaz_name(atts$basehaz))
    cat("\n assoc:          ", paste(assoc_fmt, collapse = ", "))
  }
  cat("\n algorithm:  ", tab, atts$algorithm)
  cat("\n priors:     ", tab, "see help('prior_summary')")
  if (!is.null(atts$posterior_sample_size) && atts$algorithm == "sampling")
    cat("\n sample:     ", tab, atts$posterior_sample_size, "(posterior sample size)")
  if (mvmer) {
    obs_vals <- paste0(atts$n_yobs, " (", if (jm) "Long" else "y", 1:M, ")")
    cat("\n num obs:    ", tab, paste(obs_vals, collapse = ", "))
  }
  if (jm) {
    cat("\n num subjects:   ", atts$n_subjects)
    cat(paste0("\n num events:      ", atts$n_events, " (", 
               round(100 * atts$n_events/atts$n_subjects, 1), "%)"))
  }  
  if (!is.null(atts$n_grps))
    cat("\n groups:     ", tab, 
        paste0(names(atts$n_grps), " (", unname(atts$n_grps), ")", collapse = ", "))  
  if (atts$algorithm == "sampling") {
    maxtime <- max(atts$runtime[, "total"])
    if (maxtime == 0) maxtime <- "<0.1"
    cat("\n runtime:    ", tab, maxtime, "mins")
  } 
  
  cat("\n\nEstimates:\n")
  sel <- which(colnames(x) %in% c("mcse", "n_eff", "Rhat"))
  if (!length(sel)) {
    .printfr(x, digits)
  } else {
    xtemp <- x[, -sel, drop = FALSE]
    colnames(xtemp) <- paste(" ", colnames(xtemp))
    .printfr(xtemp, digits)
    cat("\nDiagnostics:\n")
    mcse_rhat <- format(round(x[, c("mcse", "Rhat"), drop = FALSE], digits), 
                        nsmall = digits)
    n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
    print(cbind(mcse_rhat, n_eff), quote = FALSE)
    cat("\nFor each parameter, mcse is Monte Carlo standard error, ", 
        "n_eff is a crude measure of effective sample size, ", 
        "and Rhat is the potential scale reduction factor on split chains", 
        " (at convergence Rhat=1).\n", sep = '')
  }
  invisible(x)
}

# Print string showing the type of association structure
#
# @param x A stanjm object
assoc_string <- function(x) {
  sel <- grep("^which", rownames(x$assoc), invert = TRUE, value = TRUE)
  assoc <- lapply(1:ncol(x$assoc), function(m) {
    vals <- sel[which(x$assoc[sel,m] == TRUE)]     
    paste0(vals, " (Long", m, ")")
  })
  comma(unlist(assoc))
}

