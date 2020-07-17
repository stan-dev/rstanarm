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
#' based on draws from the posterior predictive distribution. By default
#' the survival probabilities are conditional on an individual's 
#' group-specific coefficients (i.e. their individual-level random
#' effects). If prediction data is provided via the \code{newdataLong}  
#' and \code{newdataEvent} arguments, then the default behaviour is to
#' sample new group-specific coefficients for the individuals in the  
#' new data using a Monte Carlo scheme that conditions on their 
#' longitudinal outcome data provided in \code{newdataLong} 
#' (sometimes referred to as "dynamic predictions", see Rizopoulos
#' (2011)). This default behaviour can be stopped by specifying 
#' \code{dynamic = FALSE}, in which case the predicted survival
#' probabilities will be marginalised over the distribution of the 
#' group-specific coefficients. This has the benefit that the user does
#' not need to provide longitudinal outcome measurements for the new 
#' individuals, however, it does mean that the predictions will incorporate
#' all the uncertainty associated with between-individual variation, since
#' the predictions aren't conditional on any observed data for the individual.
#' In addition, by default, the predicted subject-specific survival 
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
#' @templateVar stanjmArg object
#' @template args-stanjm-object
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
#'   individuals -- see the description for the \code{last_time} argument below
#'   -- however also note that when generating the survival probabilities it 
#'   is of course assumed that all individuals in \code{newdataEvent} have not 
#'   yet experienced the event (that is, any variable in \code{newdataEvent} that
#'   corresponds to the event indicator will be ignored).
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
#'     original data. If \code{edist} leads to times that are beyond
#'     the maximum event or censoring time in the original data then the 
#'     estimated survival probabilities will be truncated at that point, since
#'     the estimate for the baseline hazard is not available beyond that time.}
#' }
#' @param condition A logical specifying whether the estimated 
#'     subject-specific survival probabilities at time \code{t} should be 
#'     conditioned on survival up to a fixed time point \code{u}. The default 
#'     is for \code{condition} to be set to \code{TRUE}, unless standardised survival
#'     probabilities have been requested (by specifying \code{standardise = TRUE}),
#'     in which case \code{condition} must (and will) be set to \code{FALSE}.
#'     When conditional survival probabilities are requested, the fixed
#'     time point \code{u} will be either: (i) the value specified via the 
#'     \code{last_time} argument; or if the \code{last_time} argument is 
#'     \code{NULL} then the latest observation time for each individual 
#'     (taken to be the value in the \code{times} argument if \code{newdataEvent} 
#'     is specified, or the observed event or censoring time if \code{newdataEvent} 
#'     is \code{NULL}.
#' @param last_time A scalar, character string, or \code{NULL}. This argument 
#'     specifies the last known survival time for each individual when
#'     conditional predictions are being obtained. If 
#'     \code{newdataEvent} is provided and conditional survival predictions are being
#'     obtained, then the \code{last_time} argument can be one of the following:
#'     (i) a scalar, this will use the same last time for each individual in 
#'     \code{newdataEvent}; (ii) a character string, naming a column in 
#'     \code{newdataEvent} in which to look for the last time for each individual;
#'     (iii) \code{NULL}, in which case the default is to use the time of the latest 
#'     longitudinal observation in \code{newdataLong}. If \code{newdataEvent} is
#'     \code{NULL} then the \code{last_time} argument cannot be specified 
#'     directly; instead it will be set equal to the event or censoring time for
#'     each individual in the dataset that was used to estimate the model. 
#'     If standardised survival probabilities are requested (i.e. 
#'     \code{standardise = TRUE}) then conditional survival probabilities are
#'     not allowed and therefore the \code{last_time} argument is ignored.
#' @param ids An optional vector specifying a subset of IDs for whom the 
#'   predictions should be obtained. The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdataLong} and 
#'   \code{newdataEvent} are specified, then all individuals contained in 
#'   the new data.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predictions. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided.
#' @param times A scalar, a character string, or \code{NULL}. Specifies the  
#'   times at which the estimated survival probabilities should be calculated. 
#'   It can be either: (i) \code{NULL}, in which case it will default to the last known 
#'   survival time for each individual, as determined by the \code{last_time}
#'   argument; (ii) a scalar, specifying a time to estimate the survival probability
#'   for each of the individuals; or (iii) if \code{newdataEvent} is  
#'   provided, it can be the name of a variable in \code{newdataEvent} that 
#'   indicates the time at which the survival probabilities should be calculated  
#'   for each individual. 
#' @param standardise A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged
#'   across all individuals for whom the subject-specific predictions are 
#'   being obtained. This can be used to average over the covariate and random effects
#'   distributions of the individuals used in estimating the model, or the individuals 
#'   included in the \code{newdata} arguments. This approach of
#'   averaging across the observed distribution of the covariates is sometimes
#'   referred to as a "standardised" survival curve. If \code{standardise = TRUE}, 
#'   then the \code{times} argument must be specified and it must be constant across 
#'   individuals, that is, the survival probabilities must be calculated at the 
#'   same time for all individuals.
#' @param dynamic A logical that is only relevant if new data is provided
#'   via the \code{newdataLong} and \code{newdataEvent} arguments. If 
#'   \code{dynamic = TRUE}, then new group-specific parameters are drawn for 
#'   the individuals in the new data, conditional on their longitudinal 
#'   biomarker data contained in \code{newdataLong}. These group-specific
#'   parameters are then used to generate individual-specific survival probabilities
#'   for these individuals. These are often referred to as "dynamic predictions"
#'   in the joint modelling context, because the predictions can be updated
#'   each time additional longitudinal biomarker data is collected on the individual.
#'   On the other hand, if \code{dynamic = FALSE} then the survival probabilities
#'   will just be marginalised over the distribution of the group-specific
#'   coefficients; this will mean that the predictions will incorporate all
#'   uncertainty due to between-individual variation so there will likely be
#'   very wide credible intervals on the predicted survival probabilities.
#' @param scale A scalar, specifying how much to multiply the asymptotic 
#'   variance-covariance matrix for the random effects by, which is then
#'   used as the "width" (ie. variance-covariance matrix) of the multivariate
#'   Student-t proposal distribution in the Metropolis-Hastings algorithm. This
#'   is only relevant when \code{newdataEvent} is supplied and 
#'   \code{dynamic = TRUE}, in which case new random effects are simulated
#'   for the individuals in the new data using the Metropolis-Hastings algorithm.
#' @param draws An integer indicating the number of MCMC draws to return. 
#'   The default is to set the number of draws equal to 200, or equal to the 
#'   size of the posterior sample if that is less than 200. 
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently unused.
#'
#' @note 
#'   Note that if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdataLong} and \code{newdataEvent}. This only applies if variables  
#'   were transformed before passing the data to one of the modeling functions and  
#'   \emph{not} if transformations were specified inside the model formula.
#'    
#' @return A data frame of class \code{survfit.stanjm}. The data frame includes 
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
#' @seealso \code{\link{plot.survfit.stanjm}} for plotting the estimated survival  
#'   probabilities, \code{\link{ps_check}} for for graphical checks of the estimated 
#'   survival function, and \code{\link{posterior_traj}} for estimating the
#'   marginal or subject-specific longitudinal trajectories, and 
#'   \code{\link{plot_stack_jm}} for combining plots of the estimated subject-specific
#'   longitudinal trajectory and survival function.
#'   
#' @references 
#'   Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in 
#'   joint models for longitudinal and time-to-event data. \emph{Biometrics}
#'   \strong{67}, 819.
#'      
#' @examples
#' \donttest{
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
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,15))
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
#'   ps2 <- posterior_survfit(example_jm, ids = c(7,13,15), times = 0,
#'     extrapolate = TRUE, condition = FALSE, control = list(edist = 5))
#'     
#'   # Instead we may want to estimate subject-specific survival probabilities 
#'   # for a set of new individuals. To demonstrate this, we will simply take
#'   # the first two individuals in the estimation dataset, but pass their data
#'   # via the newdata arguments so that posterior_survfit will assume we are 
#'   # predicting survival for new individuals and draw new random effects 
#'   # under a Monte Carlo scheme (see Rizopoulos (2011)).
#'   ndL <- pbcLong[pbcLong$id %in% c(1,2),]
#'   ndE <- pbcSurv[pbcSurv$id %in% c(1,2),]
#'   ps3 <- posterior_survfit(example_jm,
#'     newdataLong = ndL, newdataEvent = ndE,
#'     last_time = "futimeYears", seed = 12345)
#'   head(ps3)
#'   # We can then compare the estimated random effects for these 
#'   # individuals based on the fitted model and the Monte Carlo scheme
#'   ranef(example_jm)$Long1$id[1:2,,drop=FALSE] # from fitted model
#'   colMeans(attr(ps3, "b_new"))                # from Monte Carlo scheme
#'   
#'   # Lastly, if we wanted to obtain "standardised" survival probabilities, 
#'   # (by averaging over the observed distribution of the fixed effect 
#'   # covariates, as well as averaging over the estimated random effects
#'   # for individuals in our estimation sample or new data) then we can 
#'   # specify 'standardise = TRUE'. We can then plot the resulting 
#'   # standardised survival curve.
#'   ps4 <- posterior_survfit(example_jm, standardise = TRUE, 
#'                            times = 0, extrapolate = TRUE)
#'   plot(ps4)
#' }
#'  
posterior_survfit <- function(object, newdataLong = NULL, newdataEvent = NULL,
                              extrapolate = TRUE, control = list(), 
                              condition = NULL, last_time = NULL, prob = 0.95, 
                              ids, times = NULL, standardise = FALSE, 
                              dynamic = TRUE, scale = 1.5,
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
  dots <- list(...)
  
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
    if (!dynamic)
      stop2("Marginalised predictions for the event outcome are ",
            "not currently implemented.")
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
    ndL <- newdatas[1:M]
    ndE <- newdatas[["Event"]]   
  }
  if (!is.null(ids)) { # user specified a subset of ids
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }  
  id_list <- factor(unique(ndE[[id_var]])) # order of ids from data, not ids arg

  # Last known survival time for each individual
  if (is.null(newdataLong)) { # user did not specify newdata
    if (!is.null(last_time))
      stop("'last_time' cannot be provided when newdata is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$eventtime[as.character(id_list)]
  } else { # user specified newdata
    if (is.null(last_time)) { # use latest longitudinal observation
      max_ytimes <- do.call("cbind", lapply(ndL, function(x) 
        tapply(x[[time_var]], x[[id_var]], FUN = max)))
      last_time <- apply(max_ytimes, 1L, max)
      # re-order last-time according to id_list
      last_time <- last_time[as.character(id_list)]
    } else if (is.character(last_time) && (length(last_time) == 1L)) {
      if (!last_time %in% colnames(ndE))
        stop("Cannot find 'last_time' column named in newdataEvent.")
      last_time <- ndE[[last_time]]      
    } else if (is.numeric(last_time) && (length(last_time) == 1L)) {
      last_time <- rep(last_time, length(id_list)) 
    } else if (is.numeric(last_time) && (length(last_time) > 1L)) {
      last_time <- last_time[as.character(id_list)]
    } else {
      stop("Bug found: could not reconcile 'last_time' argument.")
    }
    names(last_time) <- as.character(id_list)
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
    ok_control_args <- c("epoints", "edist")
    control <- get_extrapolation_control(control, ok_control_args = ok_control_args)
    endtime <- if (!is.null(control$edist)) times + control$edist else maxtime
    endtime[endtime > maxtime] <- maxtime # nothing beyond end of baseline hazard 
    time_seq <- get_time_seq(control$epoints, times, endtime, simplify = FALSE)
  } else time_seq <- list(times) # no extrapolation

  # Conditional survival times
  if (is.null(condition)) {
    condition <- !standardise
  } else if (condition && standardise) {
    stop("'condition' cannot be set to TRUE if standardised survival ",
         "probabilities are requested.")
  }
  
  # Get stanmat parameter matrix for specified number of draws
  S <- posterior_sample_size(object)
  if (is.null(draws)) 
    draws <- if (S > 200L) 200L else S 
  if (draws > S)
    stop("'draws' should be <= posterior sample size (", S, ").")
  stanmat <- as.matrix(object$stanfit)
  some_draws <- isTRUE(draws < S)
  if (some_draws) {
    samp <- sample(S, draws)
    stanmat <- stanmat[samp, , drop = FALSE]
  }

  # Draw b pars for new individuals
  if (dynamic && !is.null(newdataEvent)) {
    stanmat <- simulate_b_pars(object, stanmat = stanmat, ndL = ndL, ndE = ndE,
                               ids = id_list, times = last_time, scale = scale)
    b_new <- attr(stanmat, "b_new")
    acceptance_rate <- attr(stanmat, "acceptance_rate")
  }
  
  pars <- extract_pars(object, stanmat) # list of stanmat arrays
  
  # Matrix of surv probs at each increment of the extrapolation sequence
  # NB If no extrapolation then length(time_seq) == 1L
  surv_t <- lapply(time_seq, function(t) {  
    if (!identical(length(t), length(id_list)))
      stop("Bug found: the vector of prediction times is not the same length ",
           "as the number of individuals.")
    dat <- .pp_data_jm(object, newdataLong = ndL, newdataEvent = ndE, 
                       ids = id_list, etimes = t, long_parts = FALSE)
    surv_t <- .ll_survival(object, data = dat, pars = pars, survprob = TRUE)
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
  if (condition) {
    cond_dat <- .pp_data_jm(object, newdataLong = ndL, newdataEvent = ndE, 
                            ids = id_list, etimes = last_time, long_parts = FALSE)
    # matrix of survival probs at last_time 
    cond_surv <- .ll_survival(object, data = cond_dat, pars = pars, survprob = TRUE)
    if (is.vector(cond_surv) == 1L)
      cond_surv <- t(cond_surv)        # transform if only one individual
    cond_surv[, (last_time == 0)] <- 1 # avoids possible NaN due to numerical inaccuracies
    surv <- lapply(surv_t, function(x) { # conditional survival probs
      vec <- x / cond_surv
      vec[vec > 1] <- 1 # if t was before last_time then surv prob may be > 1
      vec
    })        
  } else surv <- surv_t
  
  # Summarise posterior draws to get median and ci
  out <- do.call("rbind", lapply(
    seq_along(surv), function(x, standardise, id_list, time_seq, prob) {
      val <- median_and_bounds(surv[[x]], prob, na.rm = TRUE)
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
  
  # temporary hack so that predictive_error can call posterior_survfit
  # with two separate conditioning times...
  fn <- tryCatch(sys.call(-1)[[1]], error = function(e) NULL)
  if (!is.null(fn) && 
      grepl("predictive_error", deparse(fn), fixed = TRUE) &&
      "last_time2" %in% names(dots)) {
    last_time2 <- ndE[[dots$last_time2]]
    cond_dat2 <- .pp_data_jm(object, newdataLong = ndL, newdataEvent = ndE, 
                         ids = id_list, etimes = last_time2, long_parts = FALSE)
    cond_surv2 <- .ll_survival(object, data = cond_dat2, pars = pars, survprob = TRUE)
    if (is.vector(cond_surv2) == 1L)
      cond_surv2 <- t(cond_surv2)        # transform if only one individual
    cond_surv2[, (last_time2 == 0)] <- 1 # avoids possible NaN due to numerical inaccuracies
    surv2 <- lapply(surv_t, function(x) { # conditional survival probs
      vec <- x / cond_surv2
      vec[vec > 1] <- 1 # if t was before last_time then surv prob may be > 1
      vec
    })
    out2 <- do.call("rbind", lapply(
      seq_along(surv2), function(x, standardise, id_list, time_seq, prob) {
        val <- median_and_bounds(surv2[[x]], prob, na.rm = TRUE)
        data.frame(IDVAR = id_list, TIMEVAR = time_seq[[x]], val$med) 
      }, standardise, id_list, time_seq, prob))
    out2 <- data.frame(out2)
    colnames(out2) <- c(id_var, time_var, "survpred_eventtime")  
    out2 <- out2[order(out2[, id_var, drop = F], out2[, time_var, drop = F]), , drop = F]
    rownames(out2) <- NULL
    out <- merge(out, out2)
  }
  
  class(out) <- c("survfit.stanjm", "data.frame")
  out <- structure(out, id_var = id_var, time_var = time_var, extrapolate = extrapolate, 
            control = control, standardise = standardise, condition = condition, 
            last_time = last_time, ids = id_list, draws = draws, seed = seed, 
            offset = offset)
  if (dynamic && !is.null(newdataEvent)) {
    out <- structure(out, b_new = b_new, acceptance_rate = acceptance_rate)
  }
  out
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
#'   \code{\link[ggplot2:geom_path]{geom_line}} and used to control features
#'   of the plotted survival function.
#'      
#' @return The plot method returns a \code{ggplot} object, also of class
#'   \code{plot.survfit.stanjm}. This object can be further customised using the
#'   \pkg{ggplot2} package. It can also be passed to the function
#'   \code{plot_stack_jm}.
#'   
#' @seealso \code{\link{posterior_survfit}}, \code{\link{plot_stack_jm}},
#'   \code{\link{posterior_traj}}, \code{\link{plot.predict.stanjm}}      
#'   
#' @examples 
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific conditional survival probabilities
#'   # for all individuals in the estimation dataset.
#'   ps1 <- posterior_survfit(example_jm, extrapolate = TRUE)
#'   
#'   # We then plot the conditional survival probabilities for
#'   # a subset of individuals
#'   plot(ps1, ids = c(7,13,15))
#'   # We can change or add attributes to the plot
#'   plot(ps1, ids = c(7,13,15), limits = "none")
#'   plot(ps1, ids = c(7,13,15), xlab = "Follow up time")
#'   plot(ps1, ids = c(7,13,15), ci_geom_args = list(fill = "red"),
#'        color = "blue", linetype = 2)
#'   plot(ps1, ids = c(7,13,15), facet_scales = "fixed")
#'   
#'   # Since the returned plot is also a ggplot object, we can
#'   # modify some of its attributes after it has been returned
#'   plot1 <- plot(ps1, ids = c(7,13,15))
#'   plot1 + 
#'     ggplot2::theme(strip.background = ggplot2::element_blank()) +
#'     ggplot2::coord_cartesian(xlim = c(0, 15)) +
#'     ggplot2::labs(title = "Some plotted survival functions")
#'     
#'   # We can also combine the plot(s) of the estimated 
#'   # subject-specific survival functions, with plot(s) 
#'   # of the estimated longitudinal trajectories for the
#'   # same individuals
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,15))
#'   pt1 <- posterior_traj(example_jm, , ids = c(7,13,15))
#'   plot_surv <- plot(ps1) 
#'   plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
#'   plot_stack_jm(plot_traj, plot_surv)
#'    
#'   # Lastly, let us plot the standardised survival function
#'   # based on all individuals in our estimation dataset
#'   ps2 <- posterior_survfit(example_jm, standardise = TRUE, times = 0,
#'                           control = list(epoints = 20))
#'   plot(ps2)   
#' }
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
  if (!standardise) x$id <- factor(x[[id_var]])
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


#' @rdname plot.survfit.stanjm
#' @export
#' @importFrom ggplot2 ggplot_build facet_wrap aes_string expand_limits
#' 
#' @description The \code{plot_stack_jm} function takes arguments containing the plots of the estimated  
#' subject-specific longitudinal trajectory (or trajectories if a multivariate  
#' joint model was estimated) and the plot of the estimated subject-specific 
#' survival function and combines them into a single figure. This is most
#' easily understood by running the \strong{Examples} below.
#' 
#' @param yplot An object of class \code{plot.predict.stanjm}, returned by a 
#'   call to the generic \code{\link[=plot.predict.stanjm]{plot}} method for 
#'   objects of class \code{predict.stanjm}. If there is more than one 
#'   longitudinal outcome, then a list of such objects can be provided.
#' @param survplot An object of class \code{plot.survfit.stanjm}, returned by a 
#'   call to the generic \code{\link[=plot.survfit.stanjm]{plot}} method for 
#'   objects of class \code{survfit.stanjm}. 
#'   
#' @return \code{plot_stack_jm} returns an object of class
#'   \code{\link[bayesplot]{bayesplot_grid}} that includes plots of the
#'   estimated subject-specific longitudinal trajectories stacked on top of the 
#'   associated subject-specific survival curve.
#'   
#' @seealso \code{\link{plot.predict.stanjm}}, \code{\link{plot.survfit.stanjm}},
#'   \code{\link{posterior_predict}}, \code{\link{posterior_survfit}}
#'    
#' @examples
#' \donttest{
#'   if (!exists("example_jm")) example(example_jm)
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,15))
#'   pt1 <- posterior_traj(example_jm, ids = c(7,13,15), extrapolate = TRUE)
#'   plot_surv <- plot(ps1) 
#'   plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
#'   plot_stack_jm(plot_traj, plot_surv)
#' }
#'  
plot_stack_jm <- function(yplot, survplot) {
  
  if (!is(yplot, "list")) yplot <- list(yplot)
  
  lapply(yplot, function(x) {
    if (!is(x, "plot.predict.stanjm"))
      stop("'yplot' should be an object of class 'plot.predict.stanjm', ",
           "or a list of such objects.", call. = FALSE)
  })
  if (!is(survplot, "plot.survfit.stanjm"))
    stop("'survplot' should be an object of class 'plot.survfit.stanjm'.",
         call. = FALSE)   
  
  y_build <- lapply(yplot, ggplot_build)
  y_layout <- lapply(y_build, function(x) x$layout$panel_layout)
  y_ids <- lapply(y_layout, function(x)
    if (!"id" %in% colnames(x)) NULL else x[["id"]])
  
  e_build <- ggplot_build(survplot)
  e_layout <- e_build$layout$panel_layout    
  e_ids <- if (!"id" %in% colnames(e_layout)) NULL else e_layout[["id"]]
  
  if (!is.null(e_ids)) {
    lapply(y_ids, function(x, e_ids) {
      if (!all(sort(x) == sort(e_ids))) {
        stop("The individuals in the 'yplot' and 'survplot' appear to differ. Please ",
             "reestimate the plots using a common 'ids' argument.", call. = FALSE)
      }
    }, e_ids = e_ids)    
  }
  
  vline <- lapply(seq_along(y_build), function(m) {
    L <- length(y_build[[m]]$data)
    dat <- y_build[[m]]$data[[L]]
    if (!"xintercept" %in% colnames(dat)) {
      found <- FALSE
    } else {
      found <- TRUE
      dat <- dat[, c("PANEL", "xintercept"), drop = FALSE] 
      if (NROW(y_layout[[m]]) > 1) {
        panel_id_map <- y_layout[[m]][, c("PANEL", "id"), drop = FALSE]
        dat <- merge(dat, panel_id_map, by = "PANEL")
      }
      dat <- dat[, grep("PANEL", colnames(dat), invert = TRUE), drop = FALSE]
      colnames(dat) <- gsub("xintercept", paste0("xintercept", m), colnames(dat), fixed = TRUE)
    }
    list(dat = dat, found = found)
  })
  vline_found <- any(sapply(vline, function(x) x$found))
  if (!vline_found)
    cat("Could not find vertical line indicating last observation time in the",
        "plot of the longitudinal trajectory; you may wish to plot the longitudinal",
        "trajectories again with 'vline = TRUE' to aid interpretation.")
  vline_dat <- lapply(vline, function(x) x$dat)
  vline_alldat <- Reduce(function(...) merge(..., all = TRUE), vline_dat)
  vline_alldat$xintercept_max <- 
    apply(vline_alldat[, grep("id", colnames(vline_alldat), invert = TRUE), drop = FALSE], 1, max) 
  
  xmax <- max(sapply(c(y_build, list(e_build)), function(i) max(i$data[[1]]$x)))
  
  if ((!is.null(e_ids)) && (length(e_ids) > 20L)) {
    stop("Unable to generate 'plot_stack_jm' for this many individuals.", call. = FALSE)      
  } else if ((!is.null(e_ids)) && (length(e_ids) > 3L)) {
    warning("'plot_stack_jm' is unlikely to be legible with more than a few individuals.",
            immediate. = TRUE, call. = FALSE)
  }

  if (!is.null(e_ids)) {
    graph_facet <- facet_wrap(~ id, scales = "free", nrow = 1) 
  } else {
    graph_facet <- NULL
  }
  
  if (vline_found) {
    graph_vline <- geom_vline(aes_string(xintercept = "xintercept_max"), 
                              vline_alldat, linetype = 2)
  } else {
    graph_vline <- NULL
  }
  
  graph_xlims <- expand_limits(x = c(0, xmax))
  
  survplot_updated <- survplot + graph_xlims + graph_facet + graph_vline
 
  yplot_updated <- lapply(yplot, function(x) x + graph_xlims + graph_facet)
  
  bayesplot::bayesplot_grid(
    plots = c(yplot_updated, list(survplot_updated)), 
    grid_args = list(ncol = 1)
  )
}


# ------------------ exported but doc kept internal

#' Generic print method for \code{survfit.stanjm} objects
#' 
#' @rdname print.survfit.stanjm
#' @method print survfit.stanjm
#' @keywords internal
#' @export
#' @param x An object of class \code{survfit.stanjm}, returned by a call to 
#'   \code{\link{posterior_survfit}}.
#' @param digits Number of digits to use for formatting the time variable and 
#'   the survival probabilities.
#' @param ... Ignored.
#' 
print.survfit.stanjm <- function(x, digits = 4, ...) {
  time_var <- attr(x, "time_var")
  x <- as.data.frame(x)
  sel <- c(time_var, "survpred", "ci_lb", "ci_ub")
  for (i in sel) 
    x[[i]] <- format(round(x[[i]], digits), nsmall = digits)
  print(x, quote = FALSE)
  invisible(x)
}

# ------------------ internal

# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"
