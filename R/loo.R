# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

#' Information criteria and cross-validation
#'
#' @description For models fit using MCMC, compute approximate leave-one-out
#'   cross-validation (LOO, LOOIC) or, less preferably, the Widely Applicable
#'   Information Criterion (WAIC) using the \pkg{\link[=loo-package]{loo}}
#'   package. (For \eqn{K}-fold cross-validation see \code{\link{kfold.stanreg}}.)
#'   Functions for  model comparison, and model weighting/averaging are also
#'   provided. 
#'   
#'   \strong{Note}: these functions are not guaranteed to work
#'   properly unless the \code{data} argument was specified when the model was
#'   fit. Also, as of \pkg{loo} version \code{2.0.0} the default number of cores
#'   is now only 1, but we recommend using as many (or close to as many) cores
#'   as possible by setting the \code{cores} argument or using
#'   \code{options(mc.cores = VALUE)} to set it for an entire session.
#'
#' @aliases loo
#' @importFrom loo loo loo.function loo.matrix is.loo
#' @export
#' @template reference-loo
#' @template reference-bayesvis
#'
#' @param x For \code{loo} and \code{waic}, a fitted model object returned by
#'   one of the rstanarm modeling functions. See \link{stanreg-objects}.
#'
#'   For the \code{loo_model_weights} method, \code{x} should be a
#'   "stanreg_list" object, which is a list of fitted model objects created by
#'   \code{\link{stanreg_list}}. \code{loo_compare} also allows \code{x} to be a
#'   single stanreg object, with the remaining objects passed via \code{...}, or
#'   a single \code{stanreg_list} object.
#'
#' @param ... For \code{loo_compare.stanreg}, \code{...} can contain objects
#'   returned by the \code{loo}, \code{\link[=kfold.stanreg]{kfold}}, or
#'   \code{waic} method (see the \strong{Examples} section, below).
#'
#'   For \code{loo_model_weights}, \code{...} should contain arguments (e.g.
#'   \code{method}) to pass to the default \code{\link[loo]{loo_model_weights}}
#'   method from the \pkg{loo} package.
#'
#' @param cores,save_psis Passed to \code{\link[loo]{loo}}.
#' @param k_threshold Threshold for flagging estimates of the Pareto shape
#'   parameters \eqn{k} estimated by \code{loo}. See the \emph{How to proceed
#'   when \code{loo} gives warnings} section, below, for details.
#'
#' @return The structure of the objects returned by \code{loo} and \code{waic}
#'   methods are documented in detail in the \strong{Value} section in
#'   \code{\link[loo]{loo}} and \code{\link[loo]{waic}} (from the \pkg{loo}
#'   package).
#'
#' @section Approximate LOO CV: The \code{loo} method for stanreg objects
#'   provides an interface to the \pkg{\link[=loo-package]{loo}} package for
#'   approximate leave-one-out cross-validation (LOO). The LOO Information
#'   Criterion (LOOIC) has the same purpose as the Akaike Information Criterion
#'   (AIC) that is used by frequentists. Both are intended to estimate the
#'   expected log predictive density (ELPD) for a new dataset. However, the AIC
#'   ignores priors and assumes that the posterior distribution is multivariate
#'   normal, whereas the functions from the \pkg{loo} package do not make this
#'   distributional assumption and integrate over uncertainty in the parameters.
#'   This only assumes that any one observation can be omitted without having a
#'   major effect on the posterior distribution, which can be judged using the
#'   diagnostic plot provided by the \code{\link[loo:pareto-k-diagnostic]{plot.loo}} method and the
#'   warnings provided by the \code{\link[loo]{print.loo}} method (see the
#'   \emph{How to Use the rstanarm Package} vignette for an example of this
#'   process).
#'
#'   \subsection{How to proceed when \code{loo} gives warnings (k_threshold)}{
#'   The \code{k_threshold} argument to the \code{loo} method for \pkg{rstanarm}
#'   models is provided as a possible remedy when the diagnostics reveal
#'   problems stemming from the posterior's sensitivity to particular
#'   observations. Warnings about Pareto \eqn{k} estimates indicate observations
#'   for which the approximation to LOO is problematic (this is described in
#'   detail in Vehtari, Gelman, and Gabry (2017) and the
#'   \pkg{\link[=loo-package]{loo}} package documentation). The
#'   \code{k_threshold} argument can be used to set the \eqn{k} value above
#'   which an observation is flagged. If \code{k_threshold} is not \code{NULL}
#'   and there are \eqn{J} observations with \eqn{k} estimates above
#'   \code{k_threshold} then when \code{loo} is called it will refit the
#'   original model \eqn{J} times, each time leaving out one of the \eqn{J}
#'   problematic observations. The pointwise contributions of these observations
#'   to the total ELPD are then computed directly and substituted for the
#'   previous estimates from these \eqn{J} observations that are stored in the
#'   object created by \code{loo}. Another option to consider is K-fold
#'   cross-validation, which is documented on a separate page (see
#'   \code{\link[=kfold.stanreg]{kfold}}).
#'
#'   \strong{Note}: in the warning messages issued by \code{loo} about large
#'   Pareto \eqn{k} estimates we recommend setting \code{k_threshold} to at
#'   least \eqn{0.7}. There is a theoretical reason, explained in Vehtari,
#'   Gelman, and Gabry (2017), for setting the threshold to the stricter value
#'   of \eqn{0.5}, but in practice they find that errors in the LOO
#'   approximation start to increase non-negligibly when \eqn{k > 0.7}.
#'   }
#'
#' @seealso
#' \itemize{
#'   \item The new \href{http://mc-stan.org/loo/articles/}{\pkg{loo} package vignettes}
#'   and various \href{http://mc-stan.org/rstanarm/articles/}{\pkg{rstanarm} vignettes}
#'   for more examples using \code{loo} and related functions with \pkg{rstanarm} models.
#'   \item \code{\link[loo]{pareto-k-diagnostic}} in the \pkg{loo} package for
#'   more on Pareto \eqn{k} diagnostics.
#'   \item \code{\link{log_lik.stanreg}} to directly access the pointwise
#'   log-likelihood matrix.
#' }
#'
#' @examples
#' \donttest{
#' fit1 <- stan_glm(mpg ~ wt, data = mtcars, refresh = 0)
#' fit2 <- stan_glm(mpg ~ wt + cyl, data = mtcars, refresh = 0)
#'
#' # (for bigger models use as many cores as possible)
#' loo1 <- loo(fit1, cores = 2)
#' print(loo1)
#' loo2 <- loo(fit2, cores = 2)
#' print(loo2)
#'
#' # when comparing models the loo objects can be passed to loo_compare
#' # as individual arguments or as a list of loo objects
#' loo_compare(loo1, loo2)
#' loo_compare(list(loo1, loo2))
#' 
#' # if the fitted model objects contain a loo object in the component "loo"
#' # then the model objects can be passed directly or as a stanreg_list
#' fit1$loo <- loo1
#' fit2$loo <- loo2
#' loo_compare(fit1, fit2)
#' 
#' # if the fitted model objects contain a loo object _and_ a waic or kfold
#' # object, then the criterion argument determines which of them the comparison
#' # is based on 
#' fit1$waic <- waic(fit1)
#' fit2$waic <- waic(fit2)
#' loo_compare(fit1, fit2, criterion = "waic")
#' 
#' # the models can also be combined into a stanreg_list object, and more 
#' # informative model names can be provided to use when printing
#' model_list <- stanreg_list(fit1, fit2, model_names = c("Fewer predictors", "More predictors"))
#' loo_compare(model_list)
#'
#' fit3 <- stan_glm(mpg ~ disp * as.factor(cyl), data = mtcars, refresh = 0)
#' loo3 <- loo(fit3, cores = 2, k_threshold = 0.7)
#' loo_compare(loo1, loo2, loo3)
#'
#' # setting detail=TRUE will also print model formulas if used with
#' # loo_compare.stanreg or loo_compare.stanreg_list
#' fit3$loo <- loo3
#' model_list <- stanreg_list(fit1, fit2, fit3)
#' loo_compare(model_list, detail=TRUE)
#'
#' # Computing model weights
#' #
#' # if the objects in model_list already have 'loo' components then those
#' # will be used. otherwise loo will be computed for each model internally
#' # (in which case the 'cores' argument may also be used and is passed to loo())
#' loo_model_weights(model_list)  # defaults to method="stacking"
#' loo_model_weights(model_list,  method = "pseudobma")
#' loo_model_weights(model_list,  method = "pseudobma", BB = FALSE)
#'
#' # you can also pass precomputed loo objects directly to loo_model_weights
#' loo_list <- list(A = loo1, B = loo2, C = loo3) # names optional (affects printing)
#' loo_model_weights(loo_list)
#' }
#'
loo.stanreg <-
  function(x,
           ...,
           cores = getOption("mc.cores", 1),
           save_psis = FALSE,
           k_threshold = NULL) {
    if (model_has_weights(x))
      recommend_exact_loo(reason = "model has weights")

    user_threshold <- !is.null(k_threshold)
    if (user_threshold) {
      validate_k_threshold(k_threshold)
    } else {
      k_threshold <- 0.7
    }

    
    if (used.sampling(x)) # chain_id to pass to loo::relative_eff
      chain_id <- chain_id_for_loo(x)
    else { # ir_idx to pass to ...
      if (exists("ir_idx",x)) {
        ir_idx <- x$ir_idx
      } else if ("diagnostics" %in% names(x$stanfit@sim) &
               "ir_idx" %in% names(x$stanfit@sim$diagnostics)) {
        ir_idx <- x$stanfit@sim$diagnostics$ir_idx
      } else {
        stop("loo not available for models fit using algorithm='", x$algorithm,
             "' and importance_resampling=FALSE.", call. = FALSE)
      }
    }

    if (is.stanjm(x)) {
      ll <- log_lik(x)
      r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id, cores = cores)
      loo_x <-
        suppressWarnings(loo.matrix(
          ll,
          r_eff = r_eff,
          cores = cores,
          save_psis = save_psis
        ))
    } else if (is.stanmvreg(x)) {
      M <- get_M(x)
      ll <- do.call("cbind", lapply(1:M, function(m) log_lik(x, m = m)))
      r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id, cores = cores)
      loo_x <-
        suppressWarnings(loo.matrix(
          ll,
          r_eff = r_eff,
          cores = cores,
          save_psis = save_psis
        ))
    } else if (is_clogit(x)) {
      ll <- log_lik.stanreg(x)
      cons <- apply(ll,MARGIN = 2, FUN = function(y) sd(y) < 1e-15)
      if (any(cons)) {
        message(
          "The following strata were dropped from the ",
          "loo calculation because log-lik is constant: ",
          paste(which(cons), collapse = ", ")
        )
        ll <- ll[,!cons, drop = FALSE]
      }
      r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id, cores = cores)
      loo_x <-
        suppressWarnings(loo.matrix(
          ll,
          r_eff = r_eff,
          cores = cores,
          save_psis = save_psis
        ))
    } else {
      args <- ll_args(x)
      llfun <- ll_fun(x)
      likfun <- function(data_i, draws) {
        exp(llfun(data_i, draws))
      }
      if (used.sampling(x)) {
        r_eff <- loo::relative_eff(
          # using function method
          x = likfun,
          chain_id = chain_id,
          data = args$data,
          draws = args$draws,
          cores = cores,
          ...
        )
      } else {
        w_ir <- as.numeric(table(ir_idx))/length(ir_idx)
        ir_uidx <- which(!duplicated(ir_idx))
        draws <- args$draws
        data <- args$data
        r_eff <- pmin(sapply(1:dim(data)[1], function(i) {lik_i <- likfun(data[i,], draws)[ir_uidx]; var(lik_i)/(sum(w_ir^2*(lik_i-mean(lik_i))^2))}),length(ir_uidx))/length(ir_idx)
      }
      loo_x <- suppressWarnings(
        loo.function(
          llfun,
          data = args$data,
          draws = args$draws,
          r_eff = r_eff,
          ...,
          cores = cores,
          save_psis = save_psis
        )
      )
    }

    bad_obs <- loo::pareto_k_ids(loo_x, k_threshold)
    n_bad <- length(bad_obs)

    out <- structure(
      loo_x,
      model_name = deparse(substitute(x)),
      discrete = is_discrete(x),
      yhash = hash_y(x),
      formula = loo_model_formula(x)
    )

    if (!length(bad_obs)) {
      if (user_threshold) {
        message(
          "All pareto_k estimates below user-specified threshold of ",
          k_threshold,
          ". \nReturning loo object."
        )
      }
      return(out)
    }

    if (!user_threshold) {
      if (n_bad > 10) {
        recommend_kfold(n_bad)
      } else {
        recommend_reloo(n_bad)
      }
      return(out)
    }

    reloo_out <- reloo(x, loo_x, obs = bad_obs)
    structure(
      reloo_out,
      model_name = attr(out, "model_name"),
      discrete = attr(out, "discrete"),
      yhash = attr(out, "yhash"),
      formula = loo_model_formula(x)
    )
  }

# WAIC
#
#' @rdname loo.stanreg
#' @aliases waic
#' @importFrom loo waic waic.function waic.matrix is.waic
#' @export
#'
waic.stanreg <- function(x, ...) {
  if (!used.sampling(x))
    STOP_sampling_only("waic")
  if (is.stanjm(x)) {
    out <- waic.matrix(log_lik(x))
  } else if (is.stanmvreg(x)) {
    M <- get_M(x)
    ll <- do.call("cbind", lapply(1:M, function(m) log_lik(x, m = m)))
    out <- waic.matrix(ll)
  } else if (is_clogit(x)) {
    out <- waic.matrix(log_lik(x))
  } else {
    args <- ll_args(x)
    out <- waic.function(ll_fun(x), data = args$data, draws = args$draws)
  }
  structure(out,
            class = c("waic", "loo"),
            model_name = deparse(substitute(x)),
            discrete = is_discrete(x),
            yhash = hash_y(x),
            formula = loo_model_formula(x))
}


#' @rdname loo.stanreg
#' @aliases loo_compare
#' @importFrom loo loo_compare
#' @export
#'
#' @param detail For \code{loo_compare.stanreg} and
#'   \code{loo_compare.stanreg_list}, if \code{TRUE} then extra information
#'   about each model (currently just the model formulas) will be printed with
#'   the output.
#' @param criterion For \code{loo_compare.stanreg} and
#'   \code{loo_compare.stanreg_list}, should the comparison be based on LOO-CV
#'   (\code{criterion="loo"}), K-fold-CV (\code{criterion="kfold"}), or WAIC
#'   (\code{criterion="waic"}). The default is LOO-CV. See the \strong{Comparing
#'   models} and \strong{Examples} sections below.
#'
#' @return \code{loo_compare} returns a matrix with class 'compare.loo'. See the
#'   \strong{Comparing models} section below for more details.
#'
#' @section Comparing models: "loo" (or "waic" or "kfold") objects can be passed
#'   to the \code{\link[loo]{loo_compare}} function in the \pkg{loo} package to
#'   perform model comparison. \pkg{rstanarm} also provides a
#'   \code{loo_compare.stanreg} method that can be used if the "loo" (or "waic"
#'   or "kfold") object has been added to the fitted model object (see the
#'   \strong{Examples} section below for how to do this). This second method
#'   allows \pkg{rstanarm} to perform some extra checks that can't be done by
#'   the \pkg{loo} package itself (e.g., verifying that all models to be
#'   compared were fit using the same outcome variable).
#'
#'   \code{loo_compare} will return a matrix with one row per model and columns
#'   containing the ELPD difference and the standard error of the difference. In
#'   the first row of the matrix will be the model with the largest ELPD
#'   (smallest LOOIC) and will contain zeros (there is no difference between
#'   this model and itself). For each of the remaining models the ELPD
#'   difference and SE are reported relative to the model with the best ELPD
#'   (the first row). See the \strong{Details} section at the
#'   \code{\link[loo]{loo_compare}} page in the \pkg{loo} package for more
#'   information.
#'
loo_compare.stanreg <-
  function(x,
           ...,
           criterion = c("loo", "kfold", "waic"),
           detail = FALSE) {
    criterion <- match.arg(criterion)
    dots <- list(...)
    fits <- c(list(x), dots)
    .loo_comparison(fits, criterion = criterion, detail = detail)
  }


#' @rdname loo.stanreg
#' @export
loo_compare.stanreg_list <-
  function(x,
           ...,
           criterion = c("loo", "kfold", "waic"),
           detail = FALSE) {
    criterion <- match.arg(criterion)
    .loo_comparison(x, criterion = criterion, detail = detail)
  }

.loo_comparison <- function(fits, criterion, detail = FALSE) {
  loos <- lapply(fits, "[[", criterion)
  if (any(sapply(loos, is.null))) {
    stop("Not all objects have a ", criterion," component.", call. = FALSE)
  }
  loos <- validate_loos(loos)
  comp <- loo::loo_compare(x = loos)
  
  if (!detail) {
    formulas <- NULL
  } else {
    formulas <- lapply(loos, attr, "formula")
    names(formulas) <- sapply(loos, attr, "model_name")
  }
  
  # Note : rows of comp are ordered by ELPD, but formulas are in same order as
  # as initial order of models when passed in by user
  structure(
    comp,
    class = c("compare_rstanarm_loos", class(comp)),
    formulas = formulas, 
    criterion = criterion
  )
}

#' @keywords internal
#' @export
#' @method print compare_rstanarm_loos
print.compare_rstanarm_loos <- function(x, ...) {
  if (is.null(attr(x, "criterion"))) {
    criterion <- NA
  } else {
    criterion <- switch(
      attr(x, "criterion"),
      "loo" = "LOO-CV",
      "kfold" = "K-fold-CV",
      "waic" = "WAIC"
    )
  }
  formulas <- attr(x, "formulas")
  if (is.null(formulas) && !is.na(criterion)) {
    cat("Model comparison based on", paste0(criterion, ":"), "\n")
  } else {
    cat("Model formulas: ")
    nms <- names(formulas)
    for (j in seq_len(NROW(x))) {
      cat("\n", paste0(nms[j], ": "),
          formula_string(formulas[[j]]))
    }
    if (!is.na(criterion)) {
      cat("\n\nModel comparison based on", paste0(criterion, ":"), "\n")
    }
  }
  
  xcopy <- x
  class(xcopy) <- "compare.loo"
  print(xcopy, ...)
  
  return(invisible(x))
}


#' @rdname loo.stanreg
#' @aliases loo_model_weights
#'
#' @importFrom loo loo_model_weights
#' @export loo_model_weights
#'
#' @export
#'
#'
#' @section Model weights: The \code{loo_model_weights} method can be used to
#'   compute model weights for a \code{"stanreg_list"} object, which is a list
#'   of fitted model objects made with \code{\link{stanreg_list}}. The end of
#'   the \strong{Examples} section has a demonstration. For details see the
#'   \code{\link[loo]{loo_model_weights}} documentation in the \pkg{loo}
#'   package.
#'
loo_model_weights.stanreg_list <-
  function(x,
           ...,
           cores = getOption("mc.cores", 1),
           k_threshold = NULL) {
    
    loos <- lapply(x, function(object) object[["loo"]])
    no_loo <- sapply(loos, is.null)
    if (!any(no_loo)) {
      loo_list <- loos
    } else if (all(no_loo)) {
      message("Computing approximate LOO-CV (models do not already have 'loo' components). ")
      loo_list <- vector(mode = "list", length = length(x))
      for (j in seq_along(x)) {
        loo_list[[j]] <-
          loo.stanreg(x[[j]], cores = cores, k_threshold = k_threshold)
      }
    } else {
      stop("Found some models with 'loo' components and some without, ", 
           "but either all or none should have 'loo' components.")
    }
    wts <- loo::loo_model_weights.default(x = loo_list, ...)
    setNames(wts, names(x))
  }

# internal ----------------------------------------------------------------
validate_k_threshold <- function(k) {
  if (!is.numeric(k) || length(k) != 1) {
    stop("'k_threshold' must be a single numeric value.",
         call. = FALSE)
  } else if (k < 0) {
    stop("'k_threshold' < 0 not allowed.",
         call. = FALSE)
  } else if (k > 1) {
    warning(
      "Setting 'k_threshold' > 1 is not recommended.",
      "\nFor details see the PSIS-LOO section in help('loo-package', 'loo').",
      call. = FALSE
    )
  }
}
recommend_kfold <- function(n) {
  warning(
    "Found ", n, " observations with a pareto_k > 0.7. ",
    "With this many problematic observations we recommend calling ",
    "'kfold' with argument 'K=10' to perform 10-fold cross-validation ",
    "rather than LOO.\n",
    call. = FALSE
  )
}
recommend_reloo <- function(n) {
  warning(
    "Found ", n, " observation(s) with a pareto_k > 0.7. ",
    "We recommend calling 'loo' again with argument 'k_threshold = 0.7' ",
    "in order to calculate the ELPD without the assumption that ",
    "these observations are negligible. ", "This will refit the model ",
    n, " times to compute the ELPDs for the problematic observations directly.\n",
    call. = FALSE
  )
}
recommend_exact_loo <- function(reason) {
  stop(
    "'loo' is not supported if ", reason, ". ",
    "If refitting the model 'nobs(x)' times is feasible, ",
    "we recommend calling 'kfold' with K equal to the ",
    "total number of observations in the data to perform exact LOO-CV.\n",
    call. = FALSE
  )
}


# Refit model leaving out specific observations
#
# @param x stanreg object
# @param loo_x the result of loo(x)
# @param obs vector of observation indexes. the model will be refit length(obs)
#   times, each time leaving out one of the observations specified in 'obs'.
# @param ... unused currently
# @param refit logical, to toggle whether refitting actually happens (only used
#   to avoid refitting in tests)
#
# @return A modified version of 'loo_x'.
# @importFrom utils capture.output
reloo <- function(x, loo_x, obs, ..., refit = TRUE) {
  if (is.stanmvreg(x))
    STOP_if_stanmvreg("reloo")
  stopifnot(!is.null(x$data), is.loo(loo_x))
  
  J <- length(obs)
  d <- kfold_and_reloo_data(x)
  lls <- vector("list", J)
  message(
    J, " problematic observation(s) found.",
    "\nModel will be refit ", J, " times."
  )
  
  if (!refit)
    return(NULL)
  
  for (j in 1:J) {
    message(
      "\nFitting model ", j, " out of ", J,
      " (leaving out observation ", obs[j], ")"
    )
    omitted <- obs[j]
    
    if (is_clogit(x)) {
      strata_id <- model.weights(model.frame(x))
      omitted <- which(strata_id == strata_id[obs[j]])
    }

    if (used.optimizing(x)) {
      fit_j_call <-
        update(
          x,
          data = d[-omitted, , drop = FALSE],
          subset = rep(TRUE, nrow(d) - length(omitted)),
          evaluate = FALSE
        )
    } else {
      fit_j_call <-
        update(
          x,
          data = d[-omitted, , drop = FALSE],
          subset = rep(TRUE, nrow(d) - length(omitted)),
          evaluate = FALSE,
          refresh = 0,
          open_progress = FALSE
        )
    }
    fit_j_call$subset <- eval(fit_j_call$subset)
    fit_j_call$data <- eval(fit_j_call$data)
    if (!is.null(getCall(x)$offset)) {
      fit_j_call$offset <- x$offset[-omitted]
    }
    capture.output(
      fit_j <- suppressWarnings(eval(fit_j_call))
    )
    
    lls[[j]] <-
      log_lik.stanreg(
        fit_j,
        newdata = d[omitted, , drop = FALSE],
        offset = x$offset[omitted],
        newx = get_x(x)[omitted, , drop = FALSE],
        newz = x$z[omitted, , drop = FALSE], # NULL other than for some stan_betareg models
        stanmat = as.matrix.stanreg(fit_j)
      )
  }
  
  # compute elpd_{loo,j} for each of the held out observations
  elpd_loo <- unlist(lapply(lls, log_mean_exp))
  
  # compute \hat{lpd}_j for each of the held out observations (using log-lik
  # matrix from full posterior, not the leave-one-out posteriors)
  ll_x <- log_lik(
    object = x,
    newdata = d[obs,, drop=FALSE],
    offset = x$offset[obs]
  )
  hat_lpd <- apply(ll_x, 2, log_mean_exp)
  
  # compute effective number of parameters
  p_loo <- hat_lpd - elpd_loo
  
  # replace parts of the loo object with these computed quantities
  sel <- c("elpd_loo", "p_loo", "looic")
  loo_x$pointwise[obs, sel] <- cbind(elpd_loo, p_loo,  -2 * elpd_loo)
  loo_x$estimates[sel, "Estimate"] <- with(loo_x, colSums(pointwise[, sel]))
  loo_x$estimates[sel, "SE"] <- with(loo_x, {
    N <- nrow(pointwise)
    sqrt(N * apply(pointwise[, sel], 2, var))
  })
  loo_x$diagnostics$pareto_k[obs] <- NA
  
  return(loo_x)
}

log_sum_exp2 <- function(a,b) {
  m <- max(a,b)
  m + log(sum(exp(c(a,b) - m)))
}

# @param x numeric vector
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# log_mean_exp (just log_sum_exp(x) - log(length(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

# Get correct data to use for kfold and reloo
#
# @param x stanreg object
# @return data frame
kfold_and_reloo_data <- function(x) {
  # either data frame or environment
  d <- x[["data"]] 
  
  sub <- getCall(x)[["subset"]]
  if (!is.null(sub)) {
    keep <- eval(substitute(sub), envir = d)
  }
  
  if (is.environment(d)) {
    # make data frame
    d <- get_all_vars(formula(x), data = d) 
  } else {
    # already a data frame
    all_vars <- all.vars(formula(x))
    if (isTRUE(x$stan_function == "stan_gamm4")) {
      # see https://github.com/stan-dev/rstanarm/issues/435
      all_vars <- c(all_vars, all.vars(getCall(x)[["random"]]))
    }
    if ("." %in% all_vars) {
      all_vars <- seq_len(ncol(d))
    }
    d <- d[, all_vars, drop=FALSE]
  }
  
  if (!is.null(sub)) {
    d <- d[keep,, drop=FALSE]
  }
  
  d <- na.omit(d)
  
  if (is_clogit(x)) {
    strata_var <- as.character(getCall(x)$strata)
    d[[strata_var]] <- model.weights(model.frame(x))
  }
  
  return(d)
}


# Calculate a SHA1 hash of y
# @param x stanreg object
# @param ... Passed to digest::sha1
#
hash_y <- function(x, ...) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("Please install the 'digest' package.")
  validate_stanreg_object(x)
  y <- get_y(x)
  attributes(y) <- NULL
  digest::sha1(x = y, ...)
}

# check if discrete or continuous
# @param object stanreg object
is_discrete <- function(object) {
  if (inherits(object, "polr"))
    return(TRUE)
  if (inherits(object, "stanmvreg")) {
    fams <- fetch(family(object), "family")
    res <- sapply(fams, function(x)
      is.binomial(x) || is.poisson(x) || is.nb(x))
    return(res)
  }
  fam <- family(object)$family
  is.binomial(fam) || is.poisson(fam) || is.nb(fam)
}

# validate objects for model comparison
validate_loos <- function(loos = list()) {
  
  if (utils::packageVersion("loo") <= "2.1.0") {
    # will be checked by loo in later versions
    yhash <- lapply(loos, attr, which = "yhash")
    yhash_check <- sapply(yhash, function(x) {
      isTRUE(all.equal(x, yhash[[1]]))
    })
    if (!all(yhash_check)) {
      warning("Not all models have the same y variable.", call. = FALSE)
    }
  }
  
  discrete <- sapply(loos, attr, which = "discrete")
  if (!all(discrete == discrete[1])) {
    stop("Discrete and continuous observation models can't be compared.",
         call. = FALSE)
  }
  
  setNames(loos, nm = lapply(loos, attr, which = "model_name"))
}


# chain_id to pass to loo::relative_eff
chain_id_for_loo <- function(object) {
  dims <- dim(object$stanfit)[1:2]
  n_iter <- dims[1]
  n_chain <- dims[2]
  rep(1:n_chain, each = n_iter)
}


# model formula to store in loo object
# @param x stanreg object
loo_model_formula <- function(x) {
  form <- try(formula(x), silent = TRUE)
  if (inherits(form, "try-error") || is.null(form)) {
    form <- "formula not found"
  }
  return(form)
}



# deprecated --------------------------------------------------------------
#' @rdname loo.stanreg
#' @param loos a list of objects produced by the \code{\link{loo}} function
#' @export 
compare_models <- function(..., loos = list(), detail = FALSE) {
  .Deprecated("loo_compare")
  
  dots <- list(...)
  if (length(dots) && length(loos)) {
    stop("'...' and 'loos' can't both be specified.", call. = FALSE)
  } else if (length(dots)) {
    loos <- dots
  } else {
    stopifnot(is.list(loos))
  }
  
  loos <- validate_loos(loos)
  comp <- loo::compare(x = loos)
  structure(
    comp,
    class = c("compare_rstanarm_loos", class(comp)),
    model_names = names(loos),
    formulas = if (!detail) NULL else lapply(loos, attr, "formula")
  )
}
