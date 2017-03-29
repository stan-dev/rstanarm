# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Sam Brilleman
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

#' Estimate subject-specific or marginal survival probabilities
#' 
#' This function allows us to generate estimated survival probabilities (either 
#' subject-specific, or by marginalising over the distribution of the random effects) 
#' based on draws from the posterior predictive distribution. In both the 
#' "subject-specifc" and "marginal" situations, the predicted survival probabilities 
#' will still be \emph{conditional} on observed values of the fixed effect covariates 
#' in the longitudinal and event submodels (ie, the predictions will be obtained 
#' using either the design matrices used in the original \code{\link{stan_jm}} model
#' call, or using the covariate values provided in the \code{newdata} argument). However, 
#' if you wish to also average over the observed distribution of the fixed effect 
#' covariates then this is possible -- however we refer to these
#' as standardised survival probabilties -- see the \code{standardise} 
#' argument below.
#' 
#' @export
#' @templateVar stanjmArg object
#' @template args-stanjm-object
#' 
#' @param newdata Optionally, a new data frame in which to look 
#'   for variables with which to predict. If omitted, the model matrices are used. 
#'   If new data is provided, then it should contain the covariate values needed
#'   for all longitudinal submodel(s) and the event submodel. There is only
#'   allowed to be one row of data for each individual in \code{newdata}, that
#'   is, time-varying covariates are not allowed in the prediction dataset. Also 
#'   note that if \code{newdata} is provided, then the \code{times} argument
#'   must also be specified. See the \strong{Details} section for further 
#'   important details that need to be considered when specifying \code{newdata}.
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
#'     conditional predictions are being obtained. Should only be specified if
#'     \code{newdata} is provided, and conditional survival predictions are being
#'     obtained. A scalar will use the same last time for each individual in 
#'     \code{newdata}. A character string will name a column in \code{newdata}
#'     in which to look for the last times. If \code{last_time} is not provided then
#'     the default is to use the value provided in the \code{times} argument.} 
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
#'   If \code{newdata} is \code{NULL}, then the 
#'   \code{times} argument is optional; if it is not provided then \code{times} 
#'   will default to the last known event or censoring time for each individual,
#'   whereas if it is provided then it must be a numeric vector of length 1, and 
#'   the survival probabilities will be calculated at the same \code{times} for 
#'   all individuals in the estimation dataset.
#'   If \code{newdata} is provided, then the \code{times} argument cannot be
#'   \code{NULL}, rather, the user must provide a numeric vector of length 1 or
#'   the name of a variable in \code{newdata}, indicating the times at which 
#'   the survival probabilities should be calculated for the individuals in 
#'   \code{newdata}. 
#' @param standardise A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged
#'   across all individuals for whom the subject-specific predictions are 
#'   being obtained. This can be used to average over the covariate distribution
#'   of the individuals used in estimating the model, or the individuals 
#'   included in \code{newdata}. This approach of
#'   averaging across the observed distribution of the covariates is sometimes
#'   referred to as a "standardised" survival curve. If \code{standardise = TRUE}, 
#'   then the \code{times} argument must be specified and it must be constant across 
#'   individuals, that is, the survival probabilities must be calculated at the 
#'   same time for all individuals.
#' @param draws An integer indicating the number of MCMC draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently unused.
#'
#' @details 
#'   If the user wishes to obtain predictions using the \code{newdata}
#'   argument then several things need to be considered: \cr 
#'   \cr
#'   First, if you wish to obtain survival probabilities for "new" individuals, 
#'   meaning those who were \strong{not} part of the dataset used to estimate
#'   the model, then you will likely want to marginalise over the distribution 
#'   of the individual-level random effects.
#'   To ensure that this happens, you must ensure that the IDs provided in the 
#'   \code{id_var} column of \code{newdata} do 
#'   \strong{not} coincide with the IDs of individuals who were used in estimating  
#'   the model. Otherwise the predictions will be obtained using draws of the random  
#'   effects for the specific individual in the estimation data with the matching ID. 
#'   (Note that in the situation where you do want to obtain predictions for a given
#'   individual who was used in the estimation data but using new values for their 
#'   covariates, for example changing their treatment code or predicting at times 
#'   other than their actual observation times, then you could do this by specifying   
#'   the relevant individuals ID in the \code{id_var} columns of \code{newdata}. \cr
#'   \cr
#'   Second, if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdata}. This only applies if variables  
#'   were transformed before passing the data to one of the modeling functions and  
#'   \emph{not} if transformations were specified inside the model formula. Also  
#'   see the \strong{Note} section in \code{\link{posterior_predict}} for a note  
#'   about using the \code{newdata} argument with binomial models.
#'    
#' @return A data frame of class \code{survfit.stanjm}. The data frame includes 
#'   columns for each of the following: 
#'   (i) the median of the posterior predictions of the estimated survival
#'   probabilities (\code{survfit});
#'   (ii) each of the lower and upper limits of the corresponding uncertainty 
#'   interval for the estimated survival probabilities (\code{ci_lb} and 
#'   \code{ci_ub});
#'   (iii) a subject identifier (\code{id_var}), unless standardised survival
#'   probabilities were estimated;
#'   (iv) the time that the estimated survival probability is calculated for 
#'   (\code{time_var}).
#'   The returned object also includes a number of additional attributes.
#' 
#' @seealso \code{\link{plot.survfit.stanjm}} for plotting the estimated survival  
#'   probabilities, \code{\link{ps_check}} for for graphical checks of the estimated 
#'   survival function, and \code{\link{posterior_traj}} for estimating the
#'   marginal or subject-specific longitudinal trajectories.
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
  
  validate_stanjm_object(object)
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
    stop("posterior_survfit cannot yet be used with shared_b or shared_coef ",
         "association structures.") 
  
  # Construct prediction data
  # ndL: dataLong to be used in predictions
  # ndE: dataEvent to be used in predictions
  if (!identical(is.null(newdataLong), is.null(newdataEvent)))
    stop("Both newdataLong and newdataEvent must be supplied together.")
  if (is.null(newdataLong)) { # user did not specify newdata
    ndL <- model.frame(object)[1:M]
    ndE <- model.frame(object)$Event
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
  newpats <- if (is.null(newdataLong)) FALSE else check_pp_ids(object, id_list)
  
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
        stop("'times' cannot be NULL if newdata is specified.")
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
  
  # Last known survival time for each individual
  if (is.null(newdataLong)) { # user did not specify newdata
    if (!is.null(control$last_time))
      stop("'last_time' cannot be provided when newdata is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$eventtime[as.character(id_list)]
  } else { # user specified newdata
    user_arg <- control$last_time # can be NULL
    if (is.null(user_arg)) {
      last_time <- times
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
    draws <- S
  if (draws > S) {
    err <- paste0("'draws' should be <= posterior sample size (", S, ").")
    stop(err)
  }
  stanmat <- as.matrix(object$stanfit)
  some_draws <- isTRUE(draws < S)
  if (some_draws) {
    samp <- sample(S, draws)
    stanmat <- stanmat[samp, , drop = FALSE]
  }
  stanmat_means <- colMeans(stanmat)
  pars_means <- extract_pars(object, means = TRUE) # posterior means
  pars <- extract_pars(object, stanmat) # array of draws

  # Draw b pars for new ids
  if (newpats) {
    if (length(object$cnms) > 1L)
      stop("posterior_survfit not yet implemented for models with more than ",
           "one grouping factor.")
    len_b <- sapply(object$glmod_stuff, function(m) length(m$cnms[[id_var]]))
    #yb <- lapply(len_b, function(x) matrix(NA, draws, length(id_list) * x))
    # Log-lik function to optimise to obtain mode of new b pars
    scale <- 1.6 # scale b_vcov to determine width of proposal distribution
    b_mode <- list()
    b_vcov <- list()
    for (i in 1:length(id_list)) {
      ndL_i <- subset_ids(object, ndL, ids = id_list[[i]])
      ndE_i <- subset_ids(object, ndE, ids = id_list[[i]])
      dat_i <- jm_data(object, ndL_i, ndE_i)
      pars_i <- pars_means
      # Get modes for new b pars, used to center proposal distribution
      val <- optim(rep(0, sum(len_b)), optim_fn, object = object, 
                   stanmat = stanmat, ndL_i = ndL_i, 
                   dat_i = dat_i, pars_i = pars_i, len_b = len_b, 
                   method = "BFGS", hessian = TRUE)
      b_mode[[i]] <- val$par
      b_vcov[[i]] <- scale * solve(val$hessian)
      # Run 1:draws iterations of the MH algorithm
      #yb_i <- matrix(NA, nrow(ebeta), sum(len_b))
      #bnew_i <- list()
      #bcurrent <- b_mode[[i]]
      #for (s in 1:nrow(ebeta)) {
      #  bnew_i[[s]] <- mh_step(b = bcurrent, delta = b_mode[[i]], sigma = b_vcov[[i]], df = 4,
      #                         yppdat, eXq, assoc_parts, stanmat[samp[[s]],])
      #  yb_i[s,] <- bcurrent <- bnew_i[[s]]
      #}
    }
  }
  
  # Matrix of surv probs at each increment of the extrapolation sequence
  # NB If no extrapolation then length(time_seq) == 1L
  surv <- lapply(time_seq, function(t) {  
    if (!identical(length(t), length(id_list)))
      stop("Bug found: the vector of prediction times is not the same length ",
           "as the number of individuals.")
    dat <- jm_data(object, newdataLong = ndL, newdataEvent = ndE, 
                   ids = id_list, etimes = t, long_parts = FALSE)
    surv_t <- ll_event(data = dat, pars = pars, basehaz = basehaz, 
                       family = family, assoc = assoc, return_ll = FALSE)
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
    cond_surv <- ll_event(data = cond_dat, pars = pars, basehaz = basehaz, 
                          family = family, assoc = assoc, return_ll = FALSE)
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
      cbind(IDVAR   = if (!standardise) id_list,
            TIMEVAR = if (!standardise) time_seq[[x]] else unique(time_seq[[x]]),
            val$med, val$lb, val$ub)
      }, standardise, id_list, time_seq, prob))
  rownames(out) <- NULL
  colnames(out) <- c(if ("IDVAR" %in% colnames(out)) id_var,
                     time_var, "survpred", "ci_lb", "ci_ub")
  out <- data.frame(out)
  if (id_var %in% colnames(out)) { # data has id column -- sort by id and time
    out <- out[order(out[, id_var, drop = F], out[, time_var, drop = F]), , drop = F]
  } else { # data does not have id column -- sort by time only
    out <- out[order(out[, time_var, drop = F]), , drop = F]
  }
  class(out) <- c("survfit.stanjm", "data.frame")
  structure(out, id_var = id_var, time_var = time_var, extrapolate = extrapolate, 
            control = control, standardise = standardise, ids = id_list, 
            draws = draws, seed = seed, offset = offset)
}

optim_fn <- function(b, object, stanmat, ndL_i, dat_i, pars_i, len_b) {
  M <- get_M(object)
  stanmat_means <- t(colMeans(stanmat))
  pars_i$b <- mapply(function(b, x) {
    names(b) <- paste0("b[", x$mod_eta$Z_names, "]")
    return(t(b))
  }, b = split(b, rep(1:M, len_b)), x = dat_i$assoc_parts, SIMPLIFY = FALSE)
  ll_long_i <- lapply(1:M, function(m) {
    args <- ll_args(object, newdata = ndL_i[[m]], m = m, 
                    stanmat = stanmat_means, user_b = pars_i$b[[m]])
    fun  <- ll_fun(object, m = m)
    return(sum(sapply(seq_len(args$N), function(j) as.vector(
      fun(i = j, data = args$data[j, , drop = FALSE], draws = args$draws)))))
  })
  ll_event_i <- ll_event(dat_i, pars_i, object$basehaz, object$family, 
                         object$assoc, one_draw = TRUE)
  ll_b_i <- mvtnorm::dmvnorm(b, mean = rep(0, length(b)), 
                             sigma = VarCorr(object)[[object$id_var]], log = TRUE)
  ll_jm_i <- ll_jm(ll_long_i, ll_event_i) + ll_b_i
  return(-ll_jm_i)  
}    

mh_step <- function(b, delta, sigma, df, yppdat, eXq, assoc_parts, stanmat) {
  
  # New proposal for b pars
  b_new <- mvtnorm::rmvt(n = 1, delta = delta, sigma = sigma, df = df)
  # Calculate density for proposal distribution
  propdens <- mvtnorm::dmvt(
    x = b, delta = delta, sigma = sigma, df = df, log = TRUE)
  propdens_new <- mvtnorm::dmvt(
    x = b_new, delta = delta, sigma = sigma, df = df, log = TRUE)
  # Calculate density for target distribution
  targdens <- ll_jm(
    x = b, delta = delta, sigma = sigma, df = df, log = TRUE)
  targdens_new <- ll_jm(
    x = b_new, delta = delta, sigma = sigma, df = df, log = TRUE)
  # MH accept/reject step
  accept_ratio <- exp(targdens_new - targdens - propdens_new + propdens)
  if (accept_ratio >= runif(1)) return(b_prop) else return(b)
}

jm_data <- function(object, newdataLong = NULL, newdataEvent = NULL, 
                    ids = NULL, etimes = NULL, long_parts = TRUE, 
                    event_parts = TRUE) {
  M <- get_M(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
  ndL <- newdatas[1:M]
  ndE <- newdatas[["Event"]]   
  ndL <- subset_ids(object, ndL, ids)
  ndE <- subset_ids(object, ndE, ids)
  id_list <- unique(newdataEvent[[id_var]])
  if (!is.null(newdataEvent) && is.null(etimes)) {
    y <- eval(formula(object, m = "Event")[[2L]], newdataEvent)
    etimes  <- unclass(y)[,"time"]
    estatus <- unclass(y)[,"status"]    
  } else if (is.null(etimes)) {
    etimes  <- object$eventtime[[as.character(id_list)]]
    estatus <- object$status[[as.character(id_list)]]
  }
  res <- nlist(M, Npat = length(id_list))
  if (long_parts && event_parts) 
    lapply(newdataLong, function(x) {
      if (!time_var %in% colnames(x)) STOP_no_var(time_var)
      mt <- tapply(x[[time_var]], factor(x[[id_var]]), max)
      if (any(mt > etimes))
        stop("There appears to be observation times in the longitudinal data that ",
             "are later than the event time specified in the 'etimes' argument.")      
    }) 
  if (long_parts) {
    ydat <- lapply(1:M, function(m) pp_data(object, newdataLong[[m]], m = m))
    yX <- fetch(ydat, "X")
    yZt <- fetch(ydat, "Zt")
    yZnames <- fetch(ydat, "Znames")
    res <- c(res, nlist(yX, yZt, yZnames))
  }
  if (event_parts) {
    qnodes <- object$quadnodes
    qq <- get_quadpoints(qnodes)
    qtimes <- unlist(lapply(qq$points,  unstandardise_quadpoints,  0, etimes))
    qwts   <- unlist(lapply(qq$weights, unstandardise_quadweights, 0, etimes))
    edat <- prepare_data_table(newdataEvent, id_var, time_var)
    edat <- rolling_merge(edat, ids = rep(id_list, qnodes), times = qtimes)
    eXq  <- .pp_data_mer_x(object, newdata = edat, m = "Event")       
    assoc_parts <- lapply(1:M, function(m) {
      ymf <- prepare_data_table(newdataLong[[m]], id_var, time_var)
      make_assoc_parts(
        ymf, assoc = object$assoc, id_var = object$id_var, 
        time_var = object$time_var, id_list = id_list, times = qtimes, 
        use_function = pp_data, object = object, m = m)
    })
    assoc_attr <- nlist(.Data = assoc_parts, qnodes, qtimes, qwts, etimes, estatus)
    assoc_parts <- do.call("structure", assoc_attr)
    res <- c(res, nlist(eXq, assoc_parts))
  }
  return(res)
}

#' Plot the estimated subject-specific or marginal survival function
#' 
#' This generic \code{plot} method for \code{survfit.stanjm} objects will
#' plot the estimated subject-specific or marginal survival function
#' using the data frame returned by a call to \code{\link{posterior_survfit}}.
#' The call to \code{posterior_survfit} should ideally have included an
#' "extrapolation" of the survival function, obtained by setting the 
#' \code{extrapolate} argument to \code{TRUE}.
#'    
#' @method plot survfit.stanjm
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
#' @param x A data frame and object of class \code{survfit.stanjm}
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
#' @return A \code{ggplot} object, also of class \code{plot.survfit.stanjm}.
#'   This object can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack}}.
#'   
#' @seealso \code{\link{posterior_survfit}}, \code{\link{plot_stack}},
#'   \code{\link{posterior_traj}}, \code{\link{plot.predict.stanjm}}      
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
plot.survfit.stanjm <- function(x, ids = NULL, 
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
      stop("The following 'ids' are not present in the survfit.stanjm object: ",
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
  class(ret) <- c("plot.survfit.stanjm", class_ret)
  ret
}


# internal ----------------------------------------------------------------

ll_args.stanjm <- function(object, newdata, m = 1, stanmat = NULL, user_b = NULL,
                           reloo_or_kfold = calling_fun %in% c("kfold", "reloo"), ...) {
  validate_stanjm_object(object)
  f <- family(object, m = m)
  draws <- nlist(f)
  has_newdata <- !is.null(newdata)
  if (model_has_weights(object))
    STOP_if_stanjm("posterior_survfit with weights")
  
  if (has_newdata) {
    ppdat <- pp_data(object, as.data.frame(newdata), offset = offset, m = m)
    x <- ppdat$x
    y <- eval(formula(object, m = m)[[2L]], newdata)
    z <- t(ppdat$Zt)
  } else {
    x <- get_x(object, m = m)
    y <- get_y(object, m = m)
    z <- get_z(object, m = m)
  }
  if (is.null(stanmat)) 
    stanmat <- as.matrix.stanreg(object)
  one_draw <- is.vector(stanmat)
  
  fname <- f$family
  if (!is.binomial(fname)) {
    data <- data.frame(y, x)
  } else {
    if (NCOL(y) == 2L) {
      trials <- rowSums(y)
      y <- y[, 1L]
    } else {
      trials <- 1
      if (is.factor(y)) 
        y <- fac2bin(y)
      stopifnot(all(y %in% c(0, 1)))
    }
    data <- data.frame(y, trials, x)
  }
  data <- cbind(data, as.matrix(z))
  
  if (one_draw) { # stanmat is a vector
    nms <- collect_nms(names(stanmat), get_M(object))
    draws$beta <- stanmat[nms$y[[m]]]
    m_stub <- get_m_stub(m)
    if (is.gaussian(fname)) 
      draws$sigma <- stanmat[[paste0(m_stub, "sigma")]]
    if (is.gamma(fname)) 
      draws$shape <- stanmat[[paste0(m_stub, "shape")]]
    if (is.ig(fname)) 
      draws$lambda <- stanmat[[paste0(m_stub, "lambda")]]
    if (is.nb(fname)) 
      draws$size <- stanmat[[paste0(m_stub, "reciprocal_dispersion")]]
    if (is.null(user_b)) { # use b pars from stanmat
      b <- stanmat[nms$y_b[[m]]]
      if (has_newdata) {
        Z_names <- ppdat$Z_names
        if (is.null(Z_names)) {
          b <- b[!grepl("_NEW_", names(b), fixed = TRUE)]
        } else {
          b <- t(as.matrix(b))
          b <- pp_b_ord(b, Z_names)
          b <- unlist(as.data.frame(b)) # keeps colnames
        }
      }
    } else { # b pars provided directly
      if (!is.vector(user_b))
        stop("'user_b' should be the same form as stanmat, i.e. a vector.")
      b <- user_b
    }
    draws$beta <- c(draws$beta, b)
  } else { # stanmat is a matrix
    nms <- collect_nms(colnames(stanmat), get_M(object))
    draws$beta <- stanmat[, nms$y[[m]], drop = FALSE]
    m_stub <- get_m_stub(m)
    if (is.gaussian(fname)) 
      draws$sigma <- stanmat[, paste0(m_stub, "sigma")]
    if (is.gamma(fname)) 
      draws$shape <- stanmat[, paste0(m_stub, "shape")]
    if (is.ig(fname)) 
      draws$lambda <- stanmat[, paste0(m_stub, "lambda")]
    if (is.nb(fname)) 
      draws$size <- stanmat[, paste0(m_stub, "reciprocal_dispersion")] 
    if (is.null(user_b)) { # use b pars from stanmat
      b <- stanmat[, nms$y_b[[m]], drop = FALSE]
      if (has_newdata) {
        Z_names <- ppdat$Z_names
        if (is.null(Z_names)) {
          b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
        } else {
          b <- pp_b_ord(b, Z_names)
        }
      }
    } else { # b pars provided directly
      if (!any(is.matrix(user_b), is.array(user_b)))
        stop("'user_b' should be the same form as stanmat, i.e. a matrix or array.")
      if (!nrow(user_b) == nrow(stanmat))
        stop("'user_b' should have the same number of rows as stanmat.")
      b <- user_b
    }
    draws$beta <- cbind(draws$beta, b)
  }
  
  nlist(data, draws, S = NROW(draws$beta), N = nrow(data))
}

# Return survival probability or log-likelihood for event submodel
#
# @param pars A list of parameter estimates, being a single draw, with elements
#   $ybeta, $ebeta, $abeta, $basehaz_coefs, $yb
# @param eXq Design matrix for event submodel, evaluated at etimes (if not NULL)
#   and qtimes
# @param basehaz A named list containing information about the baseline hazard
# @param assoc An array with information about the desired association structure
# @param assoc_parts A named list with the design matrices etc for evaluating
#   the longitudinal submodel quantities that are used in the association
#   structure
# @param family A list of family objects for each longitudinal submodel
# @param qnodes An integer specifying the number of quadrature nodes
# @param qtimes A vector of quadrature times
# @param qwts A vector of quadrature weights corresponding to the qtimes, already
#   incorporating the (b-a)/2 scaling
# @param etimes A vector of event times. If not NULL then the first length(etimes)
#   rows of eXq correspond to the etimes, and the last length(qtimes) rows 
#   correspond to the qtimes.
# @param estatus A vector of event indicators corresponding to the etimes
# @param one_draw A logical specifying whether the parameters provided in the 
#   pars argument are vectors for a single realisation of the parameter (e.g.
#   a single MCMC draw, or a posterior mean) (TRUE) or a stanmat array (FALSE)
# @param return_ll A logical specifying whether to return the log likelihood for 
#   the event submodel (TRUE) or the survival probability (FALSE)
# @param A vector with length equal to the number of individuals
ll_event <- function(data, pars, basehaz, family = NULL, assoc = NULL, 
                     one_draw = FALSE, return_ll = TRUE) {
  etimes  <- attr(data$assoc_parts, "etimes")
  estatus <- attr(data$assoc_parts, "estatus")
  qnodes  <- attr(data$assoc_parts, "qnodes")
  qtimes  <- attr(data$assoc_parts, "qtimes")
  qwts    <- attr(data$assoc_parts, "qwts")
  
  # Linear predictor for the event submodel
  e_eta <- linear_predictor(pars$ebeta, data$eXq) 
  if (one_draw) {
    aXq <- make_assoc_terms(parts = data$assoc_parts, assoc = assoc, 
                            family = family, beta = pars$beta, b = pars$b)
    e_eta <- e_eta + linear_predictor.default(pars$abeta, aXq)
  } else {
    aXq <- matrix(NA, NROW(data$eXq), NCOL(pars$abeta))
    for (s in 1:NROW(e_eta)) {
      abeta_s <- pars$abeta[s,]
      beta_s  <- lapply(pars$beta, function(x) x[s,])
      b_s     <- lapply(pars$b,    function(x) x[s,])
      aXq_s   <- make_assoc_terms(parts = data$assoc_parts, assoc = assoc, 
                                  family = family, beta = beta_s, b = b_s)
      e_eta[s,] <- e_eta[s,] + linear_predictor.default(abeta_s, aXq_s)
    }
  }

  # Baseline hazard
  if (basehaz$type_name == "weibull") { # pars$bhcoef == weibull shape
    log_basehaz <- as.vector(log(pars$bhcoef)) + 
      linear_predictor(pars$bhcoef - 1, log(qtimes))
  } else if (basehaz$type_name == "bs") { # pars$bhcoef == spline coefs
    log_basehaz <- linear_predictor(pars$bhcoef, predict(basehaz$bs_basis, qtimes))
  } else {
    stop("Not yet implemented for basehaz = ", basehaz$type_name)
  }  
  loghaz <- log_basehaz + e_eta # log haz at etimes (if not NULL) and qtimes
  
  # Calculate survival prob or log_lik  
  if (one_draw) {
    qhaz <- tail(exp(loghaz), length(qtimes)) # haz at qtimes
    qwhaz <- qwts * qhaz
    splitting_vec <- rep(1:qnodes, each = data$Npat)
    cumhaz <- Reduce('+', split(qwhaz, splitting_vec))
  } else {
    qhaz <- exp(loghaz[, tail(1:ncol(loghaz), length(qtimes))])
    qwhaz <- t(apply(qhaz, 1L, function(row) qwts * row))
    cumhaz <- Reduce('+', array2list(qwhaz, nsplits = qnodes))  
  }
  ll_survt <- -cumhaz
  if (!return_ll) { # return surv prob at time t (upper limit of integral)
    return(exp(ll_survt)) 
  } else { # return log_lik at event time
    if (is.null(etimes) || is.null(estatus))
      stop("'etimes' and 'estatus' cannot be NULL if 'return_ll = TRUE'.")
    if (one_draw) { # return vector of length npat
      return(estatus * head(loghaz, length(etimes)) + ll_survt)
    } else { # return S * npat matrix
      eloghaz <- loghaz[, 1:length(etimes), drop = FALSE]
      ll_hazt <- t(apply(eloghaz, 1L, function(row) estatus * row))
      return(ll_hazt + ll_survt)
    }
  }
} 

# Return a data.table with the key set using the appropriate time variable
# 
# @param data A data frame
# @param id_var The name of the ID variable 
# @param time_var The name of the time variable
# @return A data.table (which will be used in a rolling merge against the
#   event times and/or quadrature times)
prepare_data_table <- function(data, id_var, time_var) {
  if (survival::is.Surv(data[[1]])) { # event submodel model.frame from fitted JM
    # If the design matrix is for the event submodel and data is model.frame obtained 
    # from the fitted model (ie, not supplied by the user) then the time point
    # for merging on covariate values is taken to be either:
    # (i) the unique observation time (single row per individual surv data), or 
    # (ii) "start" of the start/stop interval (multiple row per individual surv data)
    resp_type <- attr(data[[1]], "type")
    data <- cbind(unclass(data[[1]]), data[,-1])
    if (resp_type == "right") { # single row data
      data <- data.table::data.table(data, key = c(id_var, "time"))
      data[["time"]] <- as.numeric(data[["time"]])
    } else if (resp_type == "counting") { # start/stop multiple row data
      data <- data.table::data.table(data, key = c(id_var, "start"))
      data[["start"]] <- as.numeric(data[["start"]])
    } else {
      stop("Bug found: 'data' arg appears to be the model.frame from the fitted model, ",
           "but cannot find an appropriate time variable in the Surv(.) response.")
    }
  } else { # user provided new data, or the data is for the long. submodel
    # Alternatively, the user provided the new data for the event submodel 
    # which must be single row per individual (since multiple row per individual
    # data is not allowed by posterior_survfit), or, the data is for a 
    # longitudinal submodel. In either case, the time_var variable is used
    # for merging on covariate values -- if time_var doesn't already exist in
    # the single row per individual data then we create a dud time_var variable.
    if (!time_var %in% colnames(data)) 
      data[[time_var]] <- rep(0.0, nrow(data))  
    data <- data.table::data.table(data, key = c(id_var, time_var))
  }
  return(data)
}

# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"
