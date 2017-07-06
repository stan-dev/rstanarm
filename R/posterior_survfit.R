# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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

#' Estimate subject-specific or standardised survival probabilities
#' 
#' This function allows us to generate estimated survival probabilities 
#' based on draws from the posterior predictive distribution. If 
#' \code{newdataLong} and \code{newdataEvent} are provided, then the 
#' group-specific coefficients for the new individuals
#' will be drawn using a Monte Carlo scheme conditional on the new longitudinal
#' outcome data (sometimes referred to as "dynamic predictions", see Rizopoulos
#' (2011)). In addition, by default the predicted subject-specific survival 
#' probabilities are conditional on observed values of the fixed effect 
#' covariates (ie, the predictions will be obtained using either the design 
#' matrices used in the original \code{\link{stan_jm}} model call, or using the 
#' covariate values provided in the \code{newdataLong} and \code{newdataEvent} 
#' arguments). However, if you wish to average over the observed distribution 
#' of the fixed effect covariates then this is possible -- such predictions
#' are sometimes referred to as standardised survival probabilties -- see the 
#' \code{standardise} argument below.
#' 
#' @export
#' @templateVar stanmvregArg object
#' @template args-stanmvreg-object
#' 
#' @param newdataLong,newdataEvent Optionally, a data frame (or in the case of 
#'   \code{newdataLong} this can be a list of data frames) in which to look 
#'   for variables with which to predict. If omitted, the model matrices are used. 
#'   If new data is provided, then it should also contain the longitudinal 
#'   outcome data on which to condition when drawing the new group-specific 
#'   coefficients for individuals in the new data. Note that there is only
#'   allowed to be one row of data for each individual in \code{newdataEvent}, 
#'   that is, time-varying covariates are not allowed in the prediction data for
#'   the event submodel. Also, \code{newdataEvent} can optionally include a 
#'   variable with information about the last known survival time for the new
#'   individuals -- see the description for the \code{control} argument below
#'   -- however also note that it is assumed that all individuals in
#'   \code{newdataEvent} have not yet experienced the event.
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   survival probabilities beyond the times specified in the \code{times} argument.
#'   If \code{TRUE} then the extrapolation can be further controlled using
#'   the \code{control} argument.
#' @param control A named list with parameters controlling extrapolation 
#'   of the estimated survival function when \code{extrapolate = TRUE}. The list
#'   can contain one or more of the following named elements: \cr
#'   \describe{
#'     \item{\code{epoints}}{a positive integer specifying the number of  
#'     discrete time points at which to calculate the forecasted survival 
#'     probabilities. The default is 10.}
#'     \item{\code{edist}}{a positive scalar specifying the amount of time 
#'     across which to forecast the estimated survival function, represented 
#'     in units of the time variable \code{time_var} (from fitting the model). 
#'     The default is to extrapolate between the times specified in the 
#'     \code{times} argument and the maximum event or censoring time in the 
#'     original data. If \code{ext_distance} leads to times that are beyond
#'     the maximum event or censoring time (in the original data) then the 
#'     estimated survival probabilities will be truncated at that point, since
#'     the estimate for the baseline hazard is not available beyond that time.}
#'     \item{\code{condition}}{a logical specifying whether the estimated 
#'     subject-specific survival probabilities at time \code{t} should be 
#'     conditioned on survival up to a fixed time point \code{u}. The default 
#'     is to condition on the latest observation time for each individual 
#'     (taken to be the event or censoring time if \code{newdata} is not 
#'     specified, or the value of the \code{times} argument if \code{newdata} is 
#'     specified but no \code{last_time} is provided in the \code{control} 
#'     list, or otherwise the times provided in the \code{last_time} element
#'     of the \code{control} list).}
#'     \item{\code{last_time}}{a scalar or a character string 
#'     specifying the last known survival time for each individual for whom
#'     conditional predictions are being obtained. This is only relevant if
#'     \code{newdata} is provided, and conditional survival predictions are being
#'     obtained. A scalar will use the same last time for each individual in 
#'     \code{newdataEvent}. A character string will name a column in 
#'     \code{newdataEvent} in which to look for the last times. If \code{last_time} 
#'     is not provided then the default is to use the time of the latest 
#'     longitudinal observation.} 
#'   }
#' @param ids An optional vector specifying a subset of IDs for whom the 
#'   predictions should be obtained. The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdata} is specified,
#'   then all individuals contained in \code{newdata}.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predictions. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided.
#' @param times A numeric vector of length 1 or a character string. Specifies the  
#'   times at which to obtain the estimated survival probabilities. 
#'   If \code{times} is \code{NULL}, then it will default to the last known 
#'   event or censoring time for each individual. If \code{times} is not \code{NULL} 
#'   then it must be a numeric vector of length 1, or if \code{newdataEvent} is  
#'   provided, it can be the name of a variable in \code{newdataEvent} that 
#'   indicates the time at which the survival probabilities should be calculated  
#'   for each individual. 
#' @param standardise A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged
#'   across all individuals for whom the subject-specific predictions are 
#'   being obtained. This can be used to average over the covariate distribution
#'   of the individuals used in estimating the model, or the individuals 
#'   included in the \code{newdata} arguments. This approach of
#'   averaging across the observed distribution of the covariates is sometimes
#'   referred to as a "standardised" survival curve. If \code{standardise = TRUE}, 
#'   then the \code{times} argument must be specified and it must be constant across 
#'   individuals, that is, the survival probabilities must be calculated at the 
#'   same time for all individuals.
#' @param draws An integer indicating the number of MCMC draws to return. If 
#'   the \code{newdata} arguments are \code{NULL} then the default
#'   and maximum number of draws is the size of the posterior sample. However,
#'   if \code{newdata} is provided, then the default is to set the number of 
#'   draws equal to 200 (or equal to the size of the posterior sample if that
#'   is less than 200). This ensures that the Monte Carlo algorithm for drawing 
#'   the new group-specific coefficients doesn't take an excessive amount of time. 
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently unused.
#'
#' @note 
#'   Note that if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdataLong} and \code{newdataEvent}. This only applies if variables  
#'   were transformed before passing the data to one of the modeling functions and  
#'   \emph{not} if transformations were specified inside the model formula. Also  
#'   see the \strong{Note} section in \code{\link{posterior_predict}} for a note  
#'   about using the \code{newdataLong} argument with binomial models.
#'    
#' @return A data frame of class \code{survfit.stanmvreg}. The data frame includes 
#'   columns for each of the following: 
#'   (i) the median of the posterior predictions of the estimated survival
#'   probabilities (\code{survpred});
#'   (ii) each of the lower and upper limits of the corresponding uncertainty 
#'   interval for the estimated survival probabilities (\code{ci_lb} and 
#'   \code{ci_ub});
#'   (iii) a subject identifier (\code{id_var}), unless standardised survival
#'   probabilities were estimated;
#'   (iv) the time that the estimated survival probability is calculated for 
#'   (\code{time_var}).
#'   The returned object also includes a number of additional attributes.
#' 
#' @seealso \code{\link{plot.survfit.stanmvreg}} for plotting the estimated survival  
#'   probabilities, \code{\link{ps_check}} for for graphical checks of the estimated 
#'   survival function, and \code{\link{posterior_traj}} for estimating the
#'   marginal or subject-specific longitudinal trajectories.
#'   
#' @references 
#'   Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in 
#'   joint models for longitudinal and time-to-event data. \emph{Biometrics}
#'   \strong{67}, 819.
#'      
#' @examples
#' 
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific survival probabilities for a few
#'   # selected individuals in the estimation dataset who were  
#'   # known to survive up until their censoring time. By default
#'   # the posterior_survfit function will estimate the conditional
#'   # survival probabilities, that is, conditional on having survived
#'   # until the event or censoring time, and then by default will
#'   # extrapolate the survival predictions forward from there.  
#'   head(pbcSurv[pbcSurv$status == 0,])
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,16))
#'   head(ps1)
#'   # We can plot the estimated survival probabilities using the
#'   # associated plot function
#'   plot(ps1)
#'   
#'   # If we wanted to estimate the survival probabilities for the
#'   # same three individuals as the previous example, but this time
#'   # we won't condition on them having survived up until their 
#'   # censoring time. Instead, we will estimate their probability
#'   # of having survived between 0 and 5 years given their covariates
#'   # and their estimated random effects.
#'   # The easiest way to achieve the time scale we want (ie, 0 to 5 years)
#'   # is to specify that we want the survival time estimated at time 0
#'   # and then extrapolated forward 5 years. We also specify that we
#'   # do not want to condition on their last known survival time.
#'   ps2 <- posterior_survfit(example_jm, ids = c(7,13,16), times == 0,
#'     extrapolate = TRUE, control = list(edist = 5, condition = FALSE))
#'   ps2
#'   
#'   # Instead of estimating survival probabilities for a specific individual 
#'   # in the estimation dataset, we may want to estimate the marginal 
#'   # survival probability, that is, marginalising over the individual-level
#'   # random effects. 
#'   # Here we will estimate survival between baseline and 5 years, for a 
#'   # female who received either (i) D-penicillamine or (ii) placebo. 
#'   # To do this we will need to provide the necessary values  
#'   # of the predictors via the 'newdata' argument. However, it is important
#'   # to realise that by marginalising over the random effects 
#'   # distribution we will introduce a large amount of uncertainty into
#'   # the estimated survival probabilities. This is because we have no 
#'   # longitudinal measurements for these "new" individuals and therefore do
#'   # not have any specific information with which to estimate their random
#'   # effects. As such, there is a very wide 95% uncertainty interval 
#'   # associated with the estimated survival probabilities.
#'   nd <- data.frame(id = c("new1", "new2"),
#'                    sex = c("f", "f"), 
#'                    trt = c(1, 0))
#'   ps3 <- posterior_survfit(example_jm, newdata = nd, times = 0,
#'     extrapolate = TRUE, control = list(edist = 5, condition = FALSE))
#'   ps3
#'   
#'   # We can then plot the estimated survival functions to compare
#'   # them. To do this, we use the generic plot function.
#'   plot(ps3, limits = "none")                          
#'   
#'   # Lastly, if we wanted to obtain "standardised" survival probabilities, 
#'   # (by averaging over the observed distribution of the fixed effect 
#'   # covariates, as well as averaging over the estimated random effects
#'   # for individuals in our estimation sample) then we can specify
#'   # 'standardise = TRUE'. We can then plot the resulting standardised
#'   # survival curve.
#'   ps4 <- posterior_survfit(example_jm, standardise = TRUE, 
#'                            times = 0, extrapolate = TRUE)
#'   plot(ps4)                         
#' 
#'  
posterior_survfit <- function(object, newdataLong = NULL, newdataEvent = NULL,
                              extrapolate = TRUE, control = list(), prob = 0.95, 
                              ids, times = NULL, standardise = FALSE, 
                              draws = NULL, seed = NULL, ...) {
  validate_stanmvreg_object(object)
  if (!is.jm(object)) 
    STOP_jm_only("'posterior_survfit'")
  M        <- object$n_markers
  id_var   <- object$id_var
  time_var <- object$time_var
  basehaz  <- object$basehaz
  assoc    <- object$assoc
  family   <- family(object)
  if (!is.null(seed)) 
    set.seed(seed)
  if (missing(ids)) 
    ids <- NULL
  
  # Temporary stop, until make_assoc_terms can handle it
  sel_stop <- grep("^shared", rownames(object$assoc))
  if (any(unlist(object$assoc[sel_stop,])))
    stop("'posterior_survfit' cannot yet be used with shared_b or shared_coef ",
         "association structures.") 
  
  # Construct prediction data
  # ndL: dataLong to be used in predictions
  # ndE: dataEvent to be used in predictions
  if (!identical(is.null(newdataLong), is.null(newdataEvent)))
    stop("Both newdataLong and newdataEvent must be supplied together.")
  if (is.null(newdataLong)) { # user did not specify newdata
    dats <- get_model_data(object)
    ndL <- dats[1:M]
    ndE <- dats[["Event"]]
  } else { # user specified newdata
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
    ndL <- newdatas[1:M]
    ndE <- newdatas[["Event"]]   
  }
  if (!is.null(ids)) { # user specified a subset of ids
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }  
  id_list <- unique(ndE[[id_var]]) # order of ids from data, not ids arg
  #newpats <- if (is.null(newdataLong)) FALSE else check_pp_ids(object, id_list)
  
  # Last known survival time for each individual
  if (is.null(newdataLong)) { # user did not specify newdata
    if (!is.null(control$last_time))
      stop("'last_time' cannot be provided when newdata is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$eventtime[as.character(id_list)]
  } else { # user specified newdata
    user_arg <- control$last_time # can be NULL
    if (is.null(user_arg)) { # use latest longitudinal observation
      max_ytimes <- do.call("cbind", lapply(ndL, function(x) 
        tapply(x[[time_var]], x[[id_var]], FUN = max)))
      last_time <- apply(max_ytimes, 1L, max)
      # re-order last-time according to id_list
      last_time <- last_time[as.character(id_list)]
    } else if (is.character(user_arg) && (length(user_arg) == 1L)) {
      if (!user_arg %in% colnames(ndE))
        stop("Cannot find 'last_time' column named in newdataEvent.")
      last_time <- ndE[[user_arg]]      
    } else if (is.numeric(user_arg) && (length(user_arg) == 1L)) {
      last_time <- rep(user_arg, length(id_list)) 
    } else if (is.numeric(user_arg) && (length(user_arg) > 1L)) {
      last_time <- user_arg[as.character(id_list)]
    } else {
      stop("Bug found: could not reconcile last_time argument.")
    }
  }   
  
  # Prediction times
  if (standardise) { # standardised survival probs
    times <- 
      if (is.null(times)) {
        stop("'times' cannot be NULL for obtaining standardised survival probabilities.")
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("'times' should be a numeric vector of length 1 in order to obtain ",
             "standardised survival probabilities (the subject-specific survival ",
             "probabilities will be calculated at the specified time point, and ",
             "then averaged).")      
      }    
  } else if (is.null(newdataLong)) { # subject-specific survival probs without newdata
    times <- 
      if (is.null(times)) {
        object$eventtime[as.character(id_list)]
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("If newdata is NULL then 'times' must be NULL or a single number.")     
      }
  } else { # subject-specific survival probs with newdata
    times <- 
      if (is.null(times)) {
        times <- last_time
      } else if (is.character(times) && (length(times) == 1L)) {
        if (!times %in% colnames(ndE))
          stop("Variable specified in 'times' argument could not be found in newdata.")
        tapply(ndE[[times]], ndE[[id_var]], FUN = max)
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("If newdata is specified then 'times' can only be the name of a ",
             "variable in newdata, or a single number.")      
      }
  }
  if (!identical(length(times), length(id_list)))
    stop(paste0("length of the 'times' vector should be equal to the number of individuals ",
                "for whom predictions are being obtained (", length(id_list), ")."))     
  maxtime <- max(object$eventtime)
  if (any(times > maxtime))
    stop("'times' are not allowed to be greater than the last event or censoring ",
         "time (since unable to extrapolate the baseline hazard).")
  
  # User specified extrapolation
  if (extrapolate) {
    ok_control_args <- c("epoints", "edist", "condition", "last_time")
    control <- get_extrapolation_control(control, ok_control_args = ok_control_args,
                                         standardise = standardise)
    endtime <- if (!is.null(control$edist)) times + control$edist else maxtime
    endtime[endtime > maxtime] <- maxtime # nothing beyond end of baseline hazard 
    time_seq <- get_time_seq(control$epoints, times, endtime, simplify = FALSE)
  } else time_seq <- list(times) # no extrapolation

  # Get stanmat parameter matrix for specified number of draws
  S <- posterior_sample_size(object)
  if (is.null(draws)) 
    draws <- if (!is.null(newdataEvent) && S > 200) 200 else S 
  if (draws > S)
    stop("'draws' should be <= posterior sample size (", S, ").")
  stanmat <- as.matrix(object$stanfit)
  some_draws <- isTRUE(draws < S)
  if (some_draws) {
    samp <- sample(S, draws)
    stanmat <- stanmat[samp, , drop = FALSE]
  }
  pars_means <- extract_pars(object, means = TRUE) # list of posterior means
  pars <- extract_pars(object, stanmat) # list of stanmat arrays

  # Draw b pars for new ids
  if (!is.null(newdataEvent)) {
    if (length(object$cnms) > 1L)
      stop("posterior_survfit not yet implemented for models with more than ",
           "one grouping factor.")
    sum_p <- .p(object)[[id_var]] # total num. of b pars for each individual
    # 'scale' the asymptotic vcov of posterior, such that the scaled vcov can
    # be used as the width of the proposal distribution in the MH algorithm
    scale <- 1.6 
    # Empty matrices used to collect draws for the new b pars
    b_new <- lapply(1:length(id_list), function(x) matrix(NA, nrow(stanmat), sum_p))
    for (i in 1:length(id_list)) {
      # Design matrices for individual i only
      dat_i <- jm_data(object, ndL, ndE, etimes = last_time[[i]], ids = id_list[[i]])
      # Obtain mode and var-cov matrix of posterior distribution of new b pars
      # based on asymptotic assumptions, used as center and width of proposal
      # distribution in MH algorithm
      inits <- rep(0, sum_p)
      val <- optim(inits, optim_fn, object = object, data = dat_i, 
                   pars = pars_means, method = "BFGS", hessian = TRUE)
      delta_i <- val$par                    # asymptotic mode of posterior
      Sigma_i <- scale * solve(val$hessian) # (scaled) asymptotic vcov of posterior
      b_current <- delta_i # asympotic mode used as init value for MH algorithm
      # Run MH algorithm for each individual
      for (s in 1:nrow(stanmat)) {
        pars_s <- extract_pars(object, stanmat[s, , drop = FALSE])
        b_current <- b_new[[i]][s,] <- 
          mh_step(b_old = b_current, delta = delta_i, sigma = Sigma_i, 
                  df = 4, object = object, data = dat_i, pars = pars_s)
      }
      new_nms <- unlist(sapply(dat_i$assoc_parts, function(x) x$mod_eta$Z_names))
      colnames(b_new[[i]]) <- paste0("b[", new_nms, "]")
    }
    b_new <- do.call("cbind", b_new)      # cbind new b pars for all individuals
    b_sel <- b_names(colnames(stanmat))   
    stanmat <- stanmat[, -b_sel, drop = FALSE] # drop old b pars from stanmat
    stanmat <- cbind(stanmat, b_new)           # add new b pars to stanmat
    pars <- extract_pars(object, stanmat) # reextract pars list with new b pars
  }

  # Matrix of surv probs at each increment of the extrapolation sequence
  # NB If no extrapolation then length(time_seq) == 1L
  surv <- lapply(time_seq, function(t) {  
    if (!identical(length(t), length(id_list)))
      stop("Bug found: the vector of prediction times is not the same length ",
           "as the number of individuals.")
    dat <- jm_data(object, newdataLong = ndL, newdataEvent = ndE, 
                   ids = id_list, etimes = t, long_parts = FALSE)
    surv_t <- ll_event(object, data = dat, pars = pars, survprob = TRUE)
    if (is.vector(surv_t) == 1L) 
      surv_t <- t(surv_t)   # transform if only one individual
    surv_t[, (t == 0)] <- 1 # avoids possible NaN due to numerical inaccuracies
    if (standardise) {      # standardised survival probs
      surv_t <- matrix(rowMeans(surv_t), ncol = 1)
      dimnames(surv_t) <- list(iterations = NULL, "standardised_survprob") 
    } else {
      dimnames(surv_t) <- list(iterations = NULL, ids = id_list)
    }
    surv_t
  })
  
  # If conditioning, need to obtain matrix of surv probs at last known surv time
  if (extrapolate && control$condition) {
    cond_dat <- jm_data(object, newdataLong = ndL, newdataEvent = ndE, 
                        ids = id_list, etimes = last_time, long_parts = FALSE)
    # matrix of survival probs at last_time 
    cond_surv <- ll_event(object, data = cond_dat, pars = pars, survprob = TRUE)
    if (is.vector(cond_surv) == 1L)
      cond_surv <- t(cond_surv)        # transform if only one individual
    cond_surv[, (last_time == 0)] <- 1 # avoids possible NaN due to numerical inaccuracies
    surv <- lapply(surv, function(x) { # conditional survival probs
      vec <- x / cond_surv
      vec[vec > 1] <- 1 # if t was before last_time then surv prob may be > 1
      vec
    })        
  }
  
  # Summarise posterior draws to get median and ci
  out <- do.call("rbind", lapply(
    seq_along(surv), function(x, standardise, id_list, time_seq, prob) {
      val <- median_and_bounds(surv[[x]], prob)
      if (standardise) {
        data.frame(TIMEVAR = unique(time_seq[[x]]), val$med, val$lb, val$ub)        
      } else
        data.frame(IDVAR = id_list, TIMEVAR = time_seq[[x]], val$med, val$lb, val$ub) 
      }, standardise, id_list, time_seq, prob))
  out <- data.frame(out)
  colnames(out) <- c(if ("IDVAR" %in% colnames(out)) id_var,
                     time_var, "survpred", "ci_lb", "ci_ub")  
  if (id_var %in% colnames(out)) { # data has id column -- sort by id and time
    out <- out[order(out[, id_var, drop = F], out[, time_var, drop = F]), , drop = F]
  } else { # data does not have id column -- sort by time only
    out <- out[order(out[, time_var, drop = F]), , drop = F]
  }
  rownames(out) <- NULL
  class(out) <- c("survfit.stanmvreg", "data.frame")
  structure(out, id_var = id_var, time_var = time_var, extrapolate = extrapolate, 
            control = control, standardise = standardise, ids = id_list, 
            draws = draws, seed = seed, offset = offset, 
            b_new = if (!is.null(newdataEvent)) b_new else NULL)
}

#' Plot the estimated subject-specific or marginal survival function
#' 
#' This generic \code{plot} method for \code{survfit.stanmvreg} objects will
#' plot the estimated subject-specific or marginal survival function
#' using the data frame returned by a call to \code{\link{posterior_survfit}}.
#' The call to \code{posterior_survfit} should ideally have included an
#' "extrapolation" of the survival function, obtained by setting the 
#' \code{extrapolate} argument to \code{TRUE}.
#'    
#' @method plot survfit.stanmvreg
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon 
#'   facet_wrap labs coord_cartesian
#'   
#' @templateVar idsArg ids
#' @templateVar labsArg xlab,ylab
#' @templateVar scalesArg facet_scales
#' @templateVar cigeomArg ci_geom_args
#' @template args-ids
#' @template args-labs
#' @template args-scales
#' @template args-ci-geom-args
#'  
#' @param x A data frame and object of class \code{survfit.stanmvreg}
#'   returned by a call to the function \code{\link{posterior_survfit}}.
#'   The object contains point estimates and uncertainty interval limits
#'   for estimated values of the survival function.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval for the estimated survival probability
#'   (often known as a credible interval); or \code{"none"} for no interval 
#'   limits.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_line}} and used to control features
#'   of the plotted survival function.
#'      
#' @return A \code{ggplot} object, also of class \code{plot.survfit.stanmvreg}.
#'   This object can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack}}.
#'   
#' @seealso \code{\link{posterior_survfit}}, \code{\link{plot_stack}},
#'   \code{\link{posterior_traj}}, \code{\link{plot.predict.stanmvreg}}      
#'   
#' @examples 
#' 
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific conditional survival probabilities
#'   # for all individuals in the estimation dataset.
#'   ps1 <- posterior_survfit(example_jm, extrapolate = TRUE)
#'   
#'   # We then plot the conditional survival probabilities for
#'   # a subset of individuals
#'   plot(ps1, ids = c(7,13,16))
#'   
#'   # We can change or add attributes to the plot
#'   plot(ps1, ids = c(7,13,16), limits = "none")
#'   plot(ps1, ids = c(7,13,16), xlab = "Follow up time")
#'   plot(ps1, ids = c(7,13,16), ci_geom_args = list(fill = "red"),
#'        color = "blue", linetype = 2)
#'   plot(ps1, ids = c(7,13,16), facet_scales = "fixed")
#'   
#'   # Since the returned plot is also a ggplot object, we can
#'   # modify some of its attributes after it has been returned
#'   plot1 <- plot(ps1, ids = c(7,13,16))
#'   plot1 + 
#'     ggplot2::theme(strip.background = ggplot2::element_blank()) +
#'     ggplot2::coord_cartesian(xlim = c(0, 15)) +
#'     ggplot2::labs(title = "Some plotted survival functions")
#'     
#'   # We can also combine the plot(s) of the estimated 
#'   # subject-specific survival functions, with plot(s) 
#'   # of the estimated longitudinal trajectories for the
#'   # same individuals
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,16))
#'   pt1 <- posterior_traj(example_jm, , ids = c(7,13,16))
#'   plot_surv <- plot(ps1) 
#'   plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
#'   plot_stack(plot_traj, plot_surv)
#'    
#'   # Lastly, let us plot the standardised survival function
#'   # based on all individuals in our estimation dataset
#'   ps2 <- posterior_survfit(example_jm, standardise = TRUE, times = 0,
#'                           control = list(epoints = 20))
#'   plot(ps2)   
#' 
#'    
plot.survfit.stanmvreg <- function(x, ids = NULL, 
                                limits = c("ci", "none"),  
                                xlab = NULL, ylab = NULL, facet_scales = "free", 
                                ci_geom_args = NULL, ...) {
  
  limits <- match.arg(limits)
  ci <- (limits == "ci")
  standardise <- attr(x, "standardise")
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (is.null(ylab)) ylab <- "Event free probability"
  if (!is.null(ids)) {
    if (standardise) 
      stop("'ids' argument cannot be specified when plotting standardised ",
           "survival probabilities.")
    if (!id_var %in% colnames(x))
      stop("Bug found: could not find 'id_var' column in the data frame.")
    ids_missing <- which(!ids %in% x[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the survfit.stanmvreg object: ",
           paste(ids[[ids_missing]], collapse = ", "), call. = FALSE)
    x <- x[(x[[id_var]] %in% ids), , drop = FALSE]
  } else {
    ids <- if (!standardise) attr(x, "ids") else NULL
  }
  if (!standardise) x$id <- x[[id_var]]
  x$time <- x[[time_var]]
  
  geom_defaults <- list(color = "black")
  geom_args <- set_geom_args(geom_defaults, ...)  
  
  lim_defaults <- list(alpha = 0.3)
  lim_args <- do.call("set_geom_args", c(defaults = list(lim_defaults), ci_geom_args))
  
  if ((!standardise) && (length(ids) > 60L)) {
    stop("Too many individuals to plot for. Perhaps consider limiting ",
         "the number of individuals by specifying the 'ids' argument.")
  } else if ((!standardise) && (length(ids) > 1L)) {
    graph <- ggplot(x, aes_string(x = "time", y = "survpred")) +
      theme_bw() +
      do.call("geom_line", geom_args) +
      coord_cartesian(ylim = c(0, 1)) +      
      facet_wrap(~ id, scales = facet_scales)
    if (ci) {
      lim_mapp <- list(mapping = aes_string(ymin = "ci_lb", ymax = "ci_ub"))
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  } else {
    graph <- ggplot(x, aes_string(x = "time", y = "survpred")) + 
      theme_bw() +
      do.call("geom_line", geom_args) + 
      coord_cartesian(ylim = c(0, 1))
    if (ci) {
      lim_mapp <- list(mapping = aes_string(ymin = "ci_lb", ymax = "ci_ub"))
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  }    
  
  ret <- graph + graph_limits + labs(x = xlab, y = ylab) 
  class_ret <- class(ret)
  class(ret) <- c("plot.survfit.stanmvreg", class_ret)
  ret
}


# ------------------ exported but doc kept internal

#' Generic print method for \code{survfit.stanmvreg} objects
#' 
#' @rdname print.survfit.stanmvreg
#' @method print survfit.stanmvreg
#' @keywords internal
#' @export
#' @param x An object of class \code{survfit.stanmvreg}, returned by a call to 
#'   \code{\link{posterior_survfit}}.
#' @param digits Number of digits to use for formatting the time variable and 
#'   the survival probabilities.
#' @param ... Ignored.
#' 
print.survfit.stanmvreg <- function(x, digits = 4, ...) {
  time_var <- attr(x, "time_var")
  x <- as.data.frame(x)
  sel <- c(time_var, "survpred", "ci_lb", "ci_ub")
  for (i in sel) 
    x[[i]] <- format(round(x[[i]], digits), nsmall = digits)
  print(x, quote = FALSE)
  invisible(x)
}

# ------------------ internal

# Function to optimise to obtain mode and var-cov matrix for b pars
# 
# @param b The vector of b parameters
# @param object A stanmvreg object
# @param data Output from jm_data
# @param pars Output from extract_pars
optim_fn <- function(b, object, data, pars) {
  nms <- lapply(data$assoc_parts, function(x) x$mod_eta$Z_names)
  pars <- substitute_b_pars(object, data, pars, new_b = b, new_Z_names = nms)
  ll <- ll_jm(object, data, pars, include_b = TRUE)
  return(-ll) # optimise -ll for full joint model 
}    

# Perform one iteration of the Metropolis-Hastings algorithm
# 
# @param b_old The current vector of b parameters
# @param delta The mean vector for the proposal distribution
# @param sigma The variance-covariance matrix for the proposal distribution
# @param object A stanmvreg object
# @param data Output from jm_data
# @param pars Output from extract_pars
mh_step <- function(b_old, delta, sigma, df, object, data, pars) {
  # New proposal for b vector
  b_new <- mvtnorm::rmvt(n = 1, delta = delta, sigma = sigma, df = df)
  # Calculate density for proposal distribution
  propdens_old <- mvtnorm::dmvt(x = b_old, delta, sigma, df, log = TRUE)
  propdens_new <- mvtnorm::dmvt(x = b_new, delta, sigma, df, log = TRUE)
  # Calculate density for target distribution
  nms <- lapply(data$assoc_parts, function(x) x$mod_eta$Z_names)
  pars_old <- substitute_b_pars(object, data, pars, new_b = b_old, new_Z_names = nms)
  pars_new <- substitute_b_pars(object, data, pars, new_b = b_new, new_Z_names = nms)
  targdens_old <- ll_jm(object, data, pars_old, include_b = TRUE)
  targdens_new <- ll_jm(object, data, pars_new, include_b = TRUE)
  # MH accept/reject step
  accept_ratio <- exp(targdens_new - targdens_old - propdens_new + propdens_old)
  if (accept_ratio >= runif(1)) return(b_new) else return(b_old)
}

# Function to add new b parameters
#
# @param object A stanmvreg object
# @param data Output from jm_data
# @param pars Output from extract_pars
# @param new_b A vector of new b pars, or a list of vectors with each element
#   being the new b pars for a single submodel.
# @param new_b A vector, or a list of vectors with the names for the new b pars.
substitute_b_pars <- function(object, data, pars, new_b, new_Z_names) {
  if (!is(new_b, "list")) { # split b into submodels
    len_b <- sapply(object$glmod_stuff, function(m) length(m$cnms[[object$id_var]]))
    new_b <- split(new_b, rep(1:length(len_b), len_b))
  }
  if (!is(new_Z_names, "list")) { # split Z_names into submodels
    len_b <- sapply(object$glmod_stuff, function(m) length(m$cnms[[object$id_var]]))
    new_Z_names <- split(new_Z_names, rep(1:length(len_b), len_b))
  }  
  mapply(function(x, y) {
    if (!identical(is.vector(x), is.vector(y)))
      stop("Bug found: new_b and new_Z_names should both be vectors or lists of vectors.")
    if (!identical(length(x), length(y)))
      stop("Bug found: new_b and new_Z_names should be the same length.") 
  }, new_b, new_Z_names) 
  pars$b <- mapply(function(b, nms) {
    names(b) <- paste0("b[", nms, "]")
    t(b)
  }, new_b, new_Z_names, SIMPLIFY = FALSE)
  pars$stanmat <- pars$stanmat[, -b_names(colnames(pars$stanmat)), drop = FALSE]
  pars$stanmat <- do.call("cbind", c(list(pars$stanmat), pars$b))
  return(pars)
}

# Return a data.table with the key set using the appropriate time variable
# 
# @param data A data frame
# @param id_var The name of the ID variable 
# @param time_var The name of the time variable
# @return A data.table (which will be used in a rolling merge against the
#   event times and/or quadrature times)
prepare_data_table <- function(data, id_var, time_var) {
  if (!is.data.frame(data))
    stop("'data' should be a data frame.")
  if (!id_var %in% colnames(data))
    STOP_no_var(id_var)
  if (!time_var %in% colnames(data))
    STOP_no_var(time_var)
  # ensure no rounding in data.table merge 
  data[[time_var]] <- as.numeric(data[[time_var]]) 
  data <- data.table::data.table(data, key = c(id_var, time_var))
  return(data)
}

# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"
