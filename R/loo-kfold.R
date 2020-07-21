#' K-fold cross-validation
#' 
#' The \code{kfold} method performs exact \eqn{K}-fold cross-validation. First
#' the data are randomly partitioned into \eqn{K} subsets of equal size (or as close
#' to equal as possible), or the user can specify the \code{folds} argument
#' to determine the partitioning. Then the model is refit \eqn{K} times, each time
#' leaving out one of the \eqn{K} subsets. If \eqn{K} is equal to the total
#' number of observations in the data then \eqn{K}-fold cross-validation is
#' equivalent to exact leave-one-out cross-validation (to which
#' \code{\link[=loo.stanreg]{loo}} is an efficient approximation).
#'
#' @aliases kfold
#' @importFrom loo kfold is.kfold
#' @export
#' @template reference-loo
#' 
#' @param x A fitted model object returned by one of the rstanarm modeling
#'   functions. See \link{stanreg-objects}.
#' @param K For \code{kfold}, the number of subsets (folds) into which the data
#'   will be partitioned for performing \eqn{K}-fold cross-validation. The model
#'   is refit \code{K} times, each time leaving out one of the \code{K} folds.
#'   If the \code{folds} argument is specified then \code{K} will automatically
#'   be set to \code{length(unique(folds))}, otherwise the specified value of
#'   \code{K} is passed to \code{loo::\link[loo:kfold-helpers]{kfold_split_random}} to
#'   randomly partition the data into \code{K} subsets of equal (or as close to
#'   equal as possible) size.
#' @param save_fits For \code{kfold}, if \code{TRUE}, a component \code{'fits'}
#'   is added to the returned object to store the cross-validated
#'   \link[=stanreg-objects]{stanreg} objects and the indices of the omitted
#'   observations for each fold. Defaults to \code{FALSE}.
#' @param folds For \code{kfold}, an optional integer vector with one element
#'   per observation in the data used to fit the model. Each element of the
#'   vector is an integer in \code{1:K} indicating to which of the \code{K}
#'   folds the corresponding observation belongs. There are some convenience
#'   functions available in the \pkg{loo} package that create integer vectors to
#'   use for this purpose (see the \strong{Examples} section below and also the
#'   \link[loo]{kfold-helpers} page).
#'   
#' @param cores The number of cores to use for parallelization. Instead fitting
#'   separate Markov chains for the same model on different cores, by default
#'   \code{kfold} will distribute the \code{K} models to be fit across the cores
#'   (using \code{\link[parallel:clusterApply]{parLapply}} on Windows and
#'   \code{\link[parallel]{mclapply}} otherwise). The Markov chains for each
#'   model will be run sequentially. This will often be the most efficient
#'   option, especially if many cores are available, but in some cases it may be
#'   preferable to fit the \code{K} models sequentially and instead use the
#'   cores for the Markov chains. This can be accomplished by setting
#'   \code{options(mc.cores)} to be the desired number of cores to use
#'   for the Markov chains \emph{and} also manually specifying \code{cores=1}
#'   when calling the \code{kfold} function. See the end of the
#'   \strong{Examples} section for a demonstration.
#'   
#' @param ... Currently ignored.
#'
#' @return An object with classes 'kfold' and 'loo' that has a similar structure
#'   as the objects returned by the \code{\link{loo}} and \code{\link{waic}}
#'   methods and is compatible with the \code{\link{loo_compare}} function for
#'   comparing models.
#'   
#' @examples
#' \donttest{
#' fit1 <- stan_glm(mpg ~ wt, data = mtcars, refresh = 0)
#' fit2 <- stan_glm(mpg ~ wt + cyl, data = mtcars, refresh = 0)
#' fit3 <- stan_glm(mpg ~ disp * as.factor(cyl), data = mtcars, refresh = 0)
#'
#' # 10-fold cross-validation
#' # (if possible also specify the 'cores' argument to use multiple cores)
#' (kfold1 <- kfold(fit1, K = 10))
#' kfold2 <- kfold(fit2, K = 10)
#' kfold3 <- kfold(fit3, K = 10) 
#' loo_compare(kfold1, kfold2, kfold3)
#'
#' # stratifying by a grouping variable
#' # (note: might get some divergences warnings with this model but 
#' # this is just intended as a quick example of how to code this)
#' fit4 <- stan_lmer(mpg ~ disp + (1|cyl), data = mtcars, refresh = 0)
#' table(mtcars$cyl)
#' folds_cyl <- loo::kfold_split_stratified(K = 3, x = mtcars$cyl)
#' table(cyl = mtcars$cyl, fold = folds_cyl)
#' kfold4 <- kfold(fit4, folds = folds_cyl, cores = 2)
#' print(kfold4)
#' }
#' 
#' # Example code demonstrating the different ways to specify the number 
#' # of cores and how the cores are used
#' # 
#' # options(mc.cores = NULL)
#' # 
#' # # spread the K models over N_CORES cores (method 1)
#' # kfold(fit, K, cores = N_CORES)
#' # 
#' # # spread the K models over N_CORES cores (method 2)
#' # options(mc.cores = N_CORES)
#' # kfold(fit, K)
#' #  
#' # # fit K models sequentially using N_CORES cores for the Markov chains each time
#' # options(mc.cores = N_CORES)
#' # kfold(fit, K, cores = 1)
#'
kfold.stanreg <-
  function(x,
           K = 10,
           ...,
           folds = NULL,
           save_fits = FALSE,
           cores = getOption("mc.cores", 1)) {
    
    if (is.stanmvreg(x)) {
      STOP_if_stanmvreg("kfold")
    }
    if (model_has_weights(x)) {
      stop("kfold is not currently available for models fit using weights.")
    }
    
    stopifnot(length(cores) == 1, cores == as.integer(cores), cores >= 1)
    stan_cores <- 1 
    kfold_cores <- cores
    if (kfold_cores == 1) {
      stan_cores <- getOption("mc.cores", 1)
    }
    
    
    d <- kfold_and_reloo_data(x) # defined in loo.R
    N <- nrow(d)
    
    if (is.null(folds)) {
      stopifnot(K > 1, K <= nobs(x))
      K <- as.integer(K)
      folds <- loo::kfold_split_random(K = K, N = N)
    } else {
      K <- length(unique(folds))
      stopifnot(
        length(folds) == N,
        all(folds == as.integer(folds)),
        all(folds %in% 1L:K),
        all(1:K %in% folds)
      )
      folds <- as.integer(folds)
    }
    
    calls <- list()
    omitteds <- list()
    for (k in 1:K) {
      omitted_k <- which(folds == k)
      if (used.sampling(x)) {
        fit_k_call <- update.stanreg(
            object = x,
            data = d[-omitted_k,, drop=FALSE],
            subset = rep(TRUE, nrow(d) - length(omitted_k)),
            weights = NULL,
            cores = stan_cores,
            refresh = 0,
            open_progress = FALSE,
            evaluate = FALSE # just store unevaluated calls for now
        )
      } else {
        fit_k_call <- update.stanreg(
            object = x,
            data = d[-omitted_k,, drop=FALSE],
            subset = rep(TRUE, nrow(d) - length(omitted_k)),
            weights = NULL,
            refresh = 0,
            evaluate = FALSE # just store unevaluated calls for now
        )
      }
      if (!is.null(getCall(x)$offset)) {
        fit_k_call$offset <- x$offset[-omitted_k]
      }
      fit_k_call$cores <- eval(fit_k_call$cores)
      fit_k_call$subset <- eval(fit_k_call$subset)
      fit_k_call$data <- eval(fit_k_call$data)
      fit_k_call$offset <- eval(fit_k_call$offset)
      
      omitteds[[k]] <- omitted_k
      calls[[k]] <- fit_k_call
    }
    
    
    fits <- array(list(), c(K, 2), list(NULL, c("fit", "omitted")))
    if (kfold_cores == 1) {
      lppds <- list()
      for (k in 1:K) {
        message("Fitting model ", k, " out of ", K)
        capture.output(
          fit_k <- eval(calls[[k]])
        )
        
        omitted_k <- omitteds[[k]]
        lppds[[k]] <-
          log_lik.stanreg(
            fit_k,
            newdata = d[omitted_k, , drop = FALSE],
            offset = x$offset[omitted_k],
            newx = get_x(x)[omitted_k, , drop = FALSE],
            newz = x$z[omitted_k, , drop = FALSE], # NULL other than for some stan_betareg models
            stanmat = as.matrix.stanreg(fit_k)
          )
        if (save_fits) {
          fits[k, ] <- list(fit = fit_k, omitted = omitted_k)
        }
      }
    } else { # parallelize by fold
      message("Fitting K = ", K, " models distributed over ", cores, " cores")
      if (.Platform$OS.type != "windows") {
        out <- parallel::mclapply(
          mc.cores = kfold_cores,
          mc.preschedule = FALSE,
          X = 1:K, 
          FUN = function(k) {
            fit_k <- eval(calls[[k]])
            omitted_k <- omitteds[[k]]
            lppds_k <-
              log_lik.stanreg(
                fit_k,
                newdata = d[omitted_k, , drop = FALSE],
                offset = x$offset[omitted_k],
                newx = get_x(x)[omitted_k, , drop = FALSE],
                newz = x$z[omitted_k, , drop = FALSE],
                stanmat = as.matrix.stanreg(fit_k)
              )
            return(list(lppds = lppds_k, fit = if (save_fits) fit_k else NULL))
          }
        )
      } else { # windows
        cl <- parallel::makePSOCKcluster(kfold_cores)
        on.exit(parallel::stopCluster(cl))
        out <- parallel::parLapply(
          cl = cl,
          X = 1:K,
          ...,
          fun = function(k) {
            fit_k <- eval(calls[[k]])
            omitted_k <- omitteds[[k]]
            lppds_k <-
              log_lik.stanreg(
                fit_k,
                newdata = d[omitted_k, , drop = FALSE],
                offset = x$offset[omitted_k],
                newx = get_x(x)[omitted_k, , drop = FALSE],
                newz = x$z[omitted_k, , drop = FALSE],
                stanmat = as.matrix.stanreg(fit_k)
              )
            return(list(lppds = lppds_k, fit = if (save_fits) fit_k else NULL))
          }
        )
      }
      
      lppds <- lapply(out, "[[", "lppds")
      if (save_fits) {
        for (k in 1:K) {
          fits[k, ] <- list(fit = out[[k]][["fit"]], omitted = omitteds[[k]])
        }
      }
    }
    
    elpds_unord <- unlist(lapply(lppds, function(x) {
      apply(x, 2, log_mean_exp)
    }))
    
    # make sure elpds are put back in the right order
    obs_order <- unlist(lapply(1:K, function(k) which(folds == k)))
    elpds <- rep(NA, length(elpds_unord))
    elpds[obs_order] <- elpds_unord
    
    pointwise <- cbind(elpd_kfold = elpds, p_kfold = NA, kfoldic = -2 * elpds)
    est <- colSums(pointwise)
    se_est <- sqrt(N * apply(pointwise, 2, var))
    
    out <- list(
      estimates = cbind(Estimate = est, SE = se_est),
      pointwise = pointwise,
      elpd_kfold = est[1],
      se_elpd_kfold = se_est[1]
    )
    rownames(out$estimates) <- colnames(pointwise)
    
    if (save_fits) {
      out$fits <- fits
    }

    structure(out,
              class = c("kfold", "loo"),
              K = K,
              dims = dim(lppds[[1]]),
              model_name = deparse(substitute(x)),
              discrete = is_discrete(x),
              yhash = hash_y(x),
              formula = loo_model_formula(x))
  }

