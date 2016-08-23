# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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
#
# Version Control:
# 0.0.0-1: SLB 2016-05-19:
# 0.0.0-2: SLB 2016-05-19: altered argument names in function, to be more explicit
#                          about Long and Event but they no longer align with those 
#                          in lme4, stan_glmer, etc 

#' Bayesian joint longitudinal and time-to-event models via Stan
#' 
#' Fits a shared parameter joint model for longitudinal and time-to-event 
#' (survival) data under a Bayesian framework using Stan.
#' 
#' @export
#' @templateVar armRef (Ch. 11-15)
#' @template return-stanreg-object
#' @template see-also
#' @template args-jmpriors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template reference-gelman-hill
#' 
#' @param formulaEvent A two-sided formula object describing the time-to-event
#'   submodel. The left side of the formula should be a \code{Surv()} object.
#'   See \code{\link[survival]{Surv}}
#' @param dataEvent A data frame containing the variables specified in
#'   \code{formulaEvent}.
#' @param formulaLong A two-sided linear formula object describing both the 
#'   fixed-effects and random-effects parts of the longitudinal submodel. 
#'   See \code{\link[lme4]{glmer}} for details.
#' @param dataLong A data frame containing the variables specified in
#'   \code{formulaLong}.
#' @param time_var A character string identifying the name of the variable 
#'   in \code{dataLong} which represents time.
#' @param family The family (and possibly also the link function) for the 
#'   longitudinal submodel. See \code{\link[lme4]{glmer}} for details.
#' @param base_haz A character string indicating which baseline hazard to use
#'   for the time-to-event submodel. Currently the only option allowed is 
#'   \code{"weibull"} (the default).
#' @param quadnodes A numeric scalar giving the number of quadrature nodes 
#'   for the Gauss-Kronrod quadrature which is used to approximate the 
#'   integral over the cumulative hazard in the likelihood function. Options 
#'   are 7, 11 and 15.  
#' @param subset,weights,offset Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param ... Further arguments passed to 
#'   \code{\link[rstan]{sampling}} (e.g. \code{iter}, \code{chains}, 
#'   \code{cores}, etc.) or to \code{\link[rstan]{vb}} (if \code{algorithm} is 
#'   \code{"meanfield"} or \code{"fullrank"}).
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_jm} function is used to fit a joint model (also 
#'   known as a shared parameter model) for longitudinal and time-to-event data 
#'   without needing to write the Stan model code or create the data list for Stan. 
#'   The user needs to specify the form of the longitudinal and survival submodels 
#'   using standard R formula notation as well as specify the two R data frames 
#'   which contain the data corresponding to the variables named in the formula 
#'   arguments. \cr
#'   \cr 
#'   For the longitudinal submodel a generalised linear mixed model is assumed 
#'   with any of the \code{\link[stats]{family}} choices allowed by 
#'   \code{\link[lme4]{glmer}}. For the event submodel the only option currently 
#'   available is a Weibull proportional hazards model. Time-varying covariates 
#'   are allowed in both the longitudinal submodel and the event submodel. 
#'   The association structure for the joint model can be based on any of the 
#'   following parameterisations: current expected value in the longitudinal 
#'   submodel; current value of the linear predictor in the longitudinal submodel; 
#'   first derivative (slope) in the longitudinal submodel; first derivative 
#'   (slope) for the linear predictor in the longitudinal submodel; shared random 
#'   effects; no association structure (equivalent to fitting separate longitudinal 
#'   and event models). Bayesian estimation is performed via MCMC. 
#'   The Bayesian model includes independent priors on the 
#'   regression coefficients for both the longitudinal and event submodels, 
#'   including the association parameter(s) (in much the same way as the
#'   regression parameters in \code{\link{stan_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters (in the same way as \code{\link{stan_glmer}}). 
#'   See \code{\link{priors}} for more information about the priors. \cr
#'   \cr 
#'   The \code{stan_jm.fit} function is the workhorse model fitting function 
#'   called by \code{stan_jm}.
#'   Although it is possible to call \code{stan_jm.fit} directly, it is not
#'   recommended since the user must provide several separate design matrices 
#'   in which the row orderings must correspond to the appropriate event 
#'   (or censoring) times and quadrature points. Instead, it is strongly 
#'   recommended that the user specifies the model formula and data frames 
#'   via the \code{stan_jm} function.
#' @seealso The vignette for \code{stan_glmer} and the \emph{Hierarchical 
#'   Partial Pooling} vignette.
#'    
#' @examples
#' # see help(example_model) for details on the model below
#' if (!exists("example_model")) example(example_model) 
#' print(example_model, digits = 1)
#' 
#' @import data.table
#' @importFrom survival coxph Surv
#' @importFrom lme4 glFormula glmerControl
#' @importFrom Matrix Matrix t cBind
stan_mvjm <- function(formulaEvent, dataEvent, 
                    m = 1, formulaLong, dataLong, 
                    time_var, family = list(gaussian),
                    assoc_type = list(c("currentvalue")),
                    base_haz = "weibull", quadnodes = 15, 
                    subset, weights,
                    na.action = getOption("na.action", "na.omit"),
                    offset, contrasts = NULL, m = 1, ...,
					          init = "model_based", centreLong = FALSE, centreEvent = FALSE,
					          priorLong = list(normal()), priorLong_intercept = list(normal()),
                    priorLong_ops = list(prior_options()),
                    priorEvent = normal(), priorEvent_intercept = normal(),
					          priorEvent_ops = priorEvent_options(),
					          priorAssoc = normal(),
					          priorAssoc_ops = priorAssoc_options(),
                    prior_covariance = decov(), prior_PD = FALSE, 
                    algorithm = c("sampling", "meanfield", "fullrank"),
                    adapt_delta = NULL, QR = FALSE) {


  #=============================
  # Pre-processing of arguments
  #=============================
   
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
 
  if (m < 0) stop("m must be a positive integer specifying the ",
                  "number of longitudinal markers")
      
  if ((init == "model_based") && any(unlist(c(centreLong, centreEvent)))) 
    stop("Cannot use model based initial values when 'centreLong = TRUE'",
         " or 'centreEvent = TRUE'".)
  
  if (!is.list(formulaLong)) formulaLong <- list(formulaLong)
  if (!is.list(dataLong)) dataLong <- list(dataLong)
  if (!is.list(family)) family <- list(family)
  if (!is.list(subsetLong)) subsetLong <- list(subsetLong)
  if (!is.list(weightsLong)) weightsLong <- list(weightsLong)
  if (!is.list(na.actionLong)) na.actionLong <- list(na.actionLong)
  if (!is.list(centreLong)) centreLong <- list(centreLong)
  
  unique_dataLong <- (length(dataLong) == 1)
  unique_family   <- (length(family) == 1)
  unique_subsetLong <- (length(subsetLong) == 1)
  unique_weightsLong <- (length(weightsLong) == 1)
  unique_naactionLong <- (length(na.actionLong) == 1)
  
  # Arguments that must be length m  
  arg_list <- list(formulaLong = formulaLong)                     
  sel <- (sapply(arg_list, length) == m)
  if (!all(sel)) 
    stop("The following argument(s) are lists of the incorrect length ", 
         "(should be of length ", m, "): ", 
         paste(names(arg_list)[!sel], sep = ", "))
  
  # Arguments that must be length 1 or m  
  arg_list <- list(dataLong = dataLong,
                   family = family,
                   subsetLong = subsetLong,
                   weightsLong = weightsLong,
                   na.actionLong = na.actionLong,
                   assoc_type = assoc_type,
                   priorLong = priorLong, 
                   priorLong_intercept = priorLong_intercept,
                   priorLong_ops = priorLong_ops)                     
  sel <- (sapply(arg_list, length) %in% c(1, m))
  if (!all(sel)) 
    stop("The following argument(s) are lists of the incorrect length ", 
         "(should be of length 1 or ", m, "): ", 
         paste(names(arg_list)[!sel], sep = ", "))
  
  # Set control arguments for longitudinal submodels
  controlLong <- glmerControl(check.nlev.gtreq.5 = "ignore",
                              check.nlev.gtr.1 = "stop",
                              check.nobs.vs.rankZ = "ignore",
                              check.nobs.vs.nlev = "ignore",
                              check.nobs.vs.nRE = "ignore")
  
  # Check family and link
  algorithm <- match.arg(algorithm)
  family <- lapply(family, validate_family) 
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  fam <- lapply(family, function(x) 
                which(pmatch(supported_families, x$family, nomatch = 0L) == 1L))
  if (any(lapply(fam, length) == 0L)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  
  supported_links <- lapply(fam, function(x) 
    switch(
      supported_families[x],
      binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
      gaussian = c("identity", "log", "inverse"),
      Gamma = c("identity", "log", "inverse"),
      inverse.gaussian = c("identity", "log", "inverse", "1/mu^2"),
      "neg_binomial_2" = , # intentional
      poisson = c("log", "identity", "sqrt"),
      stop("unsupported family")
    )
  )  
  link <- mapply(function(x, i) which(supported_links[[i]] == x$link),
                 family, seq_along(family), SIMPLIFY = FALSE)
  if (any(lapply(link, length) == 0L)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  #####
  if (binom_y_prop(y, family, weights))
    stop("To specify 'y' as proportion of successes and 'weights' as ",
         "number of trials please use stan_glm rather than calling ",
         "stan_glm.fit directly.")
  if (is.binomial(family$family)) {
    if (NCOL(y) == 1L) {
      if (is.numeric(y) || is.logical(y)) 
        y <- as.integer(y)
      if (is.factor(y)) 
        y <- fac2bin(y)
      if (!all(y %in% c(0L, 1L))) 
        stop("y values must be 0 or 1 for bernoulli regression.")
    } else {
      if (!isTRUE(NCOL(y) == 2L))
        stop("y should either be a vector or a matrix 1 or 2 columns.")
      trials <- as.integer(y[, 1L] + y[, 2L])
      y <- as.integer(y[, 1L])
    }
  }
  #####
                              
                       
  # Standardised GK quadrature points
  if (quadnodes == 15) {
    quadpoint.stand <- c(
      -0.991455371120812639207,
      -0.949107912342758524526,
      -0.86486442335976907279,
      -0.7415311855993944398639,
      -0.5860872354676911302941,
      -0.4058451513773971669066,
      -0.2077849550078984676007,
      0,
      0.2077849550078984676007,
      0.405845151377397166907,
      0.5860872354676911302941,
      0.741531185599394439864,
      0.86486442335976907279,
      0.9491079123427585245262,
      0.991455371120812639207) 
  } else if (quadnodes == 11) {
    quadpoint.stand <- c(
      -0.984085360094842464496,
      -0.906179845938663992798,
      -0.754166726570849220441,
      -0.5384693101056830910363,
      -0.2796304131617831934135,
      0,
      0.2796304131617831934135,
      0.5384693101056830910363,
      0.754166726570849220441,
      0.906179845938663992798,
      0.984085360094842464496)   
  } else if (quadnodes == 7) {
    quadpoint.stand <- c(
      -0.9604912687080202834235,
      -0.7745966692414833770359,
      -0.4342437493468025580021,
      0,
      0.4342437493468025580021,
      0.7745966692414833770359,
      0.9604912687080202834235)  
  } else stop("The specified number of Gauss-Kronrod quadrature points 
              ('quadnodes') must be either 7, 11 or 15.")
  
 
  #================================
  # Data for longitudinal submodel
  #================================
  
  # Items to store for each longitudinal submodel
  y             <- as.list(rep(NA,3))     # response vector
  x             <- as.list(rep(NA,3))     # design matrix with intercept
  xtemp         <- as.list(rep(NA,3))     # design matrix without intercept, possibly centred
  xbar          <- as.list(rep(NA,3))       # means of predictors
  y_centre      <- rep(0,3)        # submodel has intercept
  y_has_intercept <- rep(0,3)      # submodel has intercept
  y_N           <- rep(0,3)        # num. observations
  y_K           <- rep(0,3)        # num. predictors (excluding intercept)
  y_weights     <- as.list(rep(NA,3))     # weights
  y_offset      <- as.list(rep(NA,3))     # offsets
  Z             <- as.list(rep(NA,3))     # group related terms
  id_var        <- c()        # ID variable in each submodel
  gamma         <- c()     # initial values for intercepts
  beta          <- c()     # initial values for coefs
  
  for (i in 1:m) {
  
    if (m == 1) 
      cat("\n--> Fitting separate longitudinal model...") 
    else 
      cat(paste0("\n--> Fitting separate model for longitudinal marker ", m, "..."))  
      
    # Fit separate longitudinal model
    mod <- lme4::glmer(
             formula = formulaLong[[i]], 
             data = if (unique_dataLong) dataLong[[1]] else dataLong[[i]],
             family = if (unique_family) family[[1]] else family[[i]], 
             control = controlLong,
             subset = if (unique_subsetLong) subsetLong[[1]] else subsetLong[[i]],
             weights = if (unique_weightsLong) weightsLong[[1]] else weightsLong[[i]],
             na.action = if (unique_naactionLong) na.actionLong[[1]] else na.actionLong[[i]])
                        
    # Response vector                    
    y[[i]] <- as.vector(lme4::getME(mod, "y"))
    
    # Design matrix
    x <- as.matrix(lme4::getME(mod, "X"))
    y_has_intercept[i] <- grepl("(Intercept", colnames(x)[1L], fixed = TRUE)
    xtemp[[i]] <- if (y_has_intercept[i]) x[, -1L, drop=FALSE] else x
    
    # Centred design matrix, if required
    if (centreLong[[i]]) {
      y_centre[i] <- 1L
      xbar[[i]] <- colMeans(xtemp[[i]])
      xtemp[[i]] <- sweep(xtemp[[i]], 2, xbar[[i]], FUN = "-")
    }
    
    # Dimensions
    y_N[i] <- NROW(xtemp[[i]])
    y_K[i] <- NCOL(xtemp[[i]])
   
    # Random effect terms
    Z[[i]] <- lme4::getME(mod, "Z")
    id_var[i] <- names(mod@cnms)[1]
    y_offset[[i]] <- lme4::getME(mod, "offset")
    if (length(id_var[i]) > 1) 
      stop(paste0("Only one grouping/ID variable is currently "
                  "allowed in each longitudinal submodel."))   
    
    # Update formula if using splines or other data dependent predictors
    formulaLong[[i]] <- formula(mod)
    formvars <- grep("", attr(terms(mod), "variables"), value = TRUE)
    predvars <- grep("", attr(terms(mod), "predvars"), value = TRUE)
    if (!identical(formvars, predvars)) {
      for (j in 2:length(formvars)) {
        formulaLong[[i]] <- 
          reformulate(gsub(formvars[[j]], 
                           predvars[[j]], 
                           deparse(formulaLong[[i]][[3]]), fixed = TRUE), 
                      response = formulaLong[[i]][[2]])
      }
    }
  
    # Model based initial values
    if (init == "model_based") {
      beta[[i]] <- fixef(mod)
      if (y_has_intercept[i]) {
        gamma[[i]] <- beta[[i]][1L]
        beta[[i]] <- beta[[i]][-1L]
      } else gamma[[i]] <- 0  # not used if no intercept     
    }

  }
  sum_y_K <- sum(y_K)
  sum_y_has_intercept <- sum(y_has_intercept)
 
  # Some additional error checks
  id_var <- unique(id_var)
  if (length(id_var) != 1)
    stop("The ID variable is not the same in all longitudinal submodels.")
  id_list <- unique(lapply(group, function(x) sort(unique(x$flist))))
  if (length(id_list) != 1)
    stop("The patient IDs are not the same in all longitudinal submodels.")

  # Family indicators
  famname <- lapply(fam, function(x) supported_families[x])
  is_bernoulli  <- mapply(function(x, i)
                            is.binomial(x) && all(y[[i]] %in% 0:1),
                          famname, seq_along(famname), SIMPLIFY = FALSE)
  is_nb         <- lapply(famname, is.nb)
  is_gaussian   <- lapply(famname, is.gaussian)
  is_gamma      <- lapply(famname, is.gamma)
  is_ig         <- lapply(famname, is.ig)
  is_continuous <- lapply(seq_along(famname), function(x) 
                     (is_gaussian[[x]] || is_gamma[[x]] || is_ig[[x]]))
  
  # require intercept for certain family and link combinations
  lapply(seq_along(y_has_intercept), function(x) {
    if (!y_has_intercept[[x]]) {
      linkname <- supported_links[[x]][link[[x]]]
      needs_intercept <- 
        !is_gaussian[[x]] && linkname == "identity" ||
        is_gamma[[x]] && linkname == "inverse" ||
        is.binomial(famname[[x]]) && linkname == "log"
      if (needs_intercept)
        stop(paste0("To use the combination of family and link ", 
                    "specified for longitudinal marker ", x,
                    ", the model must have an intercept."))
    }
  })  
    
  #=========================
  # Data for event submodel
  #=========================

  # Survival submodel
  cat("\n--> Fitting separate survival model...\n") 
  
  # Set up model frame for event submodel 
  cluster_term <- paste0("cluster(", id_var, ")")
  formulaEvent <- do.call(update.formula, list(
                            formulaEvent, 
                            paste0(" ~ . +", cluster_term)))
                            
  e_mod <- survival::coxph(formula = formulaEvent, 
                                 data = dataEvent,
                                 weights = weightsEvent,
                                 subset = subsetEvent,
                                 na.action = na.actionEvent,
                                 control = controlEvent, x = TRUE)
                             
  e_mf <- model.frame(e_mod)
  e_mf <- cbind(e_mf[,1][,1:ncol(e_mf[,1])], e_mf)
  e_y <- coxmod$y
  
  # For each individual, identify final event time and event indicator
  if (attr(e_y, "type") == "counting") {
    tvc         <- TRUE
    mf_event    <- do.call(rbind, lapply(
                             split(e_mf, e_mf[, cluster_term]),
                             function(d) d[which.max(d[,"stop"]), ]))
    flist_event <- mf_event[, cluster_term]
    eventtime   <- mf_event$stop
    d           <- mf_event$status
  
    e_mf           <- data.table(cbind(e_y, e_mf), key = c(cluster_term, "start"))
    e_mf_eventtime <- e_mf[, .SD[.N], by = e_mf[, cluster_term]]
    # Unstandardised quadrature points
    quadpoint <- lapply(quadpoint.stand, FUN = function(x) 
                          (eventtime/2) * x + (eventtime/2))
    # Model frame corresponding to observation times which are 
    #   as close as possible to the unstandardised quadrature points                      
    e_mf_quadtime  <- do.call(rbind, lapply(quadpoint, FUN = function(x) 
                         e_mf[data.table::SJ(flist_event, x), 
                         roll = TRUE, rollends = c(TRUE, TRUE)]))
    # Model frame evaluated at both event times and quadrature points
    e_mf_quadtime <- rbind(e_mf_eventtime, e_mf_quadtime, idcol = "xbind.id")
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are time varying covariates in the event submodel and
    #   therefore the design matrix differs depending on the quadrature point 
    e_x_quadtime   <- update(coxmod, data = e_mf_quadtime)$x
  } else if (attr(e_y, "type") == "right") {
    tvc         <- FALSE 
    mf_event    <- e_mf
    flist_event <- mf_event[, cluster_term]
    eventtime   <- mf_event$time
    d           <- mf_event$status
    # Unstandardised quadrature points
    quadpoint <- lapply(quadpoint.stand, FUN = function(x) 
                          (eventtime/2) * x + (eventtime/2))    
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are no time varying covariates in the event submodel and
    #   therefore the design matrix is identical at event time and at all
    #   quadrature points
    e_x_quadtime   <- do.call(rbind, lapply(1:(quadnodes + 1), 
                                            FUN = function(x) coxmod$x))
  } else stop("Only 'right' or 'counting' type Surv objects are allowed 
               on the LHS of the event submodel formula")

  # Incorporate intercept term (since Cox model does not have intercept)
  e_x_quadtime  <- cbind("(Intercept)" = rep(1, NROW(e_x_quadtime)), e_x_quadtime)

  # Centering of design matrix for event model
  e_x <- as.matrix(e_x_quadtime)  
  e_has_intercept <- grepl("(Intercept", colnames(e_x)[1L], fixed = TRUE)
  e_xtemp <- if (e_has_intercept) e_x[, -1L, drop=FALSE] else e_x
  if (centreEvent) {
    e_xbar <- colMeans(e_xtemp)
    e_xtemp <- sweep(e_xtemp, 2, e_xbar, FUN = "-")
  }
  e_K <- NCOL(e_xtemp)
  Npat <- length(eventtime)
  quadweight_rep <- rep(quadweight, each = Npat)  
  eventtime_rep <- rep(eventtime, times = quadnodes)  
  quadweight_times_half_eventtime <- 0.5 * quadweight_rep * eventtime_rep   
  
  # Model based initial values
  if (init == "model_based") {
    e_beta <- c(0, e_mod$coef)
  }
    
  # Error checks for the ID variable
  if (!identical(id_list, sort(unique(flist_event))))
    stop("The patient IDs are not the same in the longitudinal and event"
         "submodels.")
 
  #====================================================================
  # Longitudinal submodel: calculate design matrices, and id vector 
  # at the event times and quadrature points
  #====================================================================
    
  # Items to store for each longitudinal submodel
  xqtemp          <- as.list(rep(NA,3))   # design matrix (without intercept) for 
                                          # longitudinal submodel calculated at event 
                                          # and quad times, possibly centred
  dxdt_quadtime   <- as.list(rep(NA,3))   # first derivative of design matrix
  Zq              <- as.list(rep(NA,3))

  # Set up a second longitudinal model frame which includes the time variable
  for (i = 1/m) {
    formulaLong_wtime <- do.call(update.formula, list(
                                    formulaLong[[i]], 
                                    paste0("~ . +", time_var)))
    mod_wtime <- lme4::glFormula(formula = formulaLong_wtime
                          data = if (unique_dataLong) dataLong[[1]] else dataLong[[i]],
                          family = if (unique_family) family[[1]] else family[[i]], 
                          control = controlLong,
                          subset = if (unique_subsetLong) subsetLong[[1]] else subsetLong[[i]],
                          weights = if (unique_weightsLong) weightsLong[[1]] else weightsLong[[i]],
                          na.action = if (unique_naactionLong) na.actionLong[[1]] else na.actionLong[[i]])
                        
    mf <- data.table(mod_wtime$fr, key = c(id_var, time_var))
  
    # Identify which row in longitudinal data is closest to event time
    mf_eventtime <- mf[data.table::SJ(flist_event, eventtime), 
                          roll = TRUE, rollends = c(FALSE, TRUE)]
  
    # Identify which row in longitudinal data is closest to quadrature point
    #   NB if the quadrature point is earlier than the first observation time, 
    #   then covariates values are carried back to avoid missing values - I
    #   should add a warning for when this is the case! In any other case, the 
    #   observed covariates values from the most recent observation time
    #   preceeding the quadrature point are carried forward to represent the 
    #   covariate value(s) at the quadrature point. (To avoid missingness  
    #   there is no limit on how far forwards or how far backwards covariate 
    #   values can be carried). If no time varying covariates are present in
    #   the longitudinal submodel (other than the time variable) then nothing 
    #   is carried forward or backward.
    mf_quadtime <- do.call(rbind, lapply(quadpoint, FUN = function(x) 
                           mf[data.table::SJ(flist_event, x), 
                           roll = TRUE, rollends = c(TRUE, TRUE)]))
  
    # Obtain long design matrix evaluated at event times (xbind.id == 1) and  
    #   quadrature points (xbind.id == 2)
    names(mf_eventtime)[names(mf_eventtime) == "eventtime"] <- time_var
    names(mf_quadtime)[names(mf_quadtime) == "quadpoint"]   <- time_var
    mf_quadtime <- rbind(mf_eventtime, mf_quadtime, idcol = "xbind.id")  
    mod_quadtime <- lme4::glFormula(formula = formulaLong
                          data = mf_quadtime,
                          family = if (unique_family) family[[1]] else family[[i]], 
                          control = controlLong)
         
    xq <- as.matrix(mod_quadtime$X)
    xqtemp[[i]] <- if (y_has_intercept[i]) xq[, -1L, drop=FALSE] else xq  
    Zq[[i]] <- getME(mod_quadtime, "Z")
      #Needs working out to appropriately deal with offsets??
      #offset_quadtime <- model.offset(mod_quadtime$fr) %ORifNULL% double(0)

    # Centering of design matrix for longitudinal model at event times
    # and quadrature times 
    if (centreLong[i]) xqtemp[[i]] <- sweep(xqtemp[[i]], 2, xbar[[i]], FUN = "-")
    
    if (any(c("currentslope", "etaslope") %in% assoc_type)) {
      #need to contruct derivative of design matrix
    } else dxdt_quadtime[[i]] <- NULL
   
  }
   
  #================================
  # Data for association structure
  #================================
  
  if (is.null(assoc_type)) {
    assoc <- 0L
    a_K <- 0L
  } else {
    assoc <- 1L
    a_K <- length(unlist(assoc_type))
  }
  
   
  #=====================
  # Prior distributions
  #=====================
 
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus")
  ok_intercept_dists <- ok_dists[1:3]

  # Priors for longitudinal submodel(s)
  priorLong_scaled <- priorLong_ops$scaled
  priorLong_min_prior_scale <- priorLong_ops$min_prior_scale
  priorLong_scale_for_dispersion <- priorLong_ops$prior_scale_for_dispersion
  
  if (is.null(priorLong)) {
    priorLong_dist <- 0L
    priorLong_mean <- as.array(rep(0, sum_y_K))
    priorLong_scale <- priorLong_df <- as.array(rep(1, sum_y_K))
  } else {
    if (!is.list(priorLong)) 
      stop("'prior' should be a named list.")
    priorLong_dist <- priorLong$dist
    priorLong_scale <- priorLong$scale
    priorLong_mean <- priorLong$location
    priorLong_df <- priorLong$df
    priorLong_df[is.na(priorLong_df)] <- 1
    if (!priorLong_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "))
    } else if (priorLong_dist %in% c("normal", "t")) {
      priorLong_dist <- ifelse(priorLong_dist == "normal", 1L, 2L)
      priorLong_scale <- set_priorLong_scale(priorLong_scale, default = 2.5, 
                                     link = family$link)
    } else {
      priorLong_dist <- ifelse(priorLong_dist == "hs", 3L, 4L)
    }
    
    priorLong_df <- maybe_broadcast(priorLong_df, sum_y_K)
    priorLong_df <- as.array(pmin(.Machine$double.xmax, priorLong_df))
    priorLong_mean <- maybe_broadcast(priorLong_mean, sum_y_K)
    priorLong_mean <- as.array(priorLong_mean)
    priorLong_scale <- maybe_broadcast(priorLong_scale, sum_y_K)
  }
  if (is.null(priorLong_intercept)) {
    priorLong_dist_for_intercept <- 0L
    priorLong_mean_for_intercept <- 0 
    priorLong_scale_for_intercept <- priorLong_df_for_intercept <- 1
  } else {
    if (!is.list(priorLong_intercept)) 
      stop("'priorLong_intercept' should be a named list.")
    priorLong_dist_for_intercept <- priorLong_intercept$dist
    priorLong_scale_for_intercept <- priorLong_intercept$scale
    priorLong_mean_for_intercept <- priorLong_intercept$location
    priorLong_df_for_intercept <- priorLong_intercept$df 
    priorLong_df_for_intercept[is.na(priorLong_df_for_intercept)] <- 1
    
    if (!priorLong_dist_for_intercept %in% unlist(ok_intercept_dists))
      stop("The prior distribution for the intercept should be one of ",
           paste(names(ok_intercept_dists), collapse = ", "))
    priorLong_dist_for_intercept <- 
      ifelse(priorLong_dist_for_intercept == "normal", 1L, 2L)
    priorLong_scale_for_intercept <- 
      set_priorLong_scale(priorLong_scale_for_intercept, default = 10, 
                      link = family$link)
    priorLong_df_for_intercept <- min(.Machine$double.xmax, priorLong_df_for_intercept)
  }

  # Priors for event submodel
  priorEvent_scaled <- priorEvent_ops$scaled
  priorEvent_min_prior_scale <- priorEvent_ops$min_prior_scale
  priorEvent_scale_for_weibull <- priorEvent_ops$prior_scale_for_weibull
  
  if (is.null(priorEvent)) {
    priorEvent_dist <- 0L
    priorEvent_mean <- as.array(rep(0, e_K))
    priorEvent_scale <- priorEvent_df <- as.array(rep(1, e_K))
  } else {
    if (!is.list(priorEvent)) 
      stop("'priorEvent' should be a named list.")
    priorEvent_dist <- priorEvent$dist
    priorEvent_scale <- priorEvent$scale
    priorEvent_mean <- priorEvent$location
    priorEvent_df <- priorEvent$df
    priorEvent_df[is.na(priorEvent_df)] <- 1
    if (!priorEvent_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the event model coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "))
    } else if (priorEvent_dist %in% c("normal", "t")) {
      priorEvent_dist <- ifelse(priorEvent_dist == "normal", 1L, 2L)
      priorEvent_scale <- set_prior_scale(priorEvent_scale, default = 2.5, 
                                     link = family$link)
    } else {
      priorEvent_dist <- ifelse(priorEvent_dist == "hs", 3L, 4L)
    }
    
    priorEvent_df <- maybe_broadcast(priorEvent_df, e_K)
    priorEvent_df <- as.array(pmin(.Machine$double.xmax, priorEvent_df))
    priorEvent_mean <- maybe_broadcast(priorEvent_mean, e_K)
    priorEvent_mean <- as.array(priorEvent_mean)
    priorEvent_scale <- maybe_broadcast(priorEvent_scale, e_K)
  }
  if (is.null(priorEvent_intercept)) {
    priorEvent_dist_for_intercept <- 0L
    priorEvent_mean_for_intercept <- 0 
    priorEvent_scale_for_intercept <- priorEvent_df_for_intercept <- 1
  } else {
    if (!is.list(priorEvent_intercept)) 
      stop("'priorEvent_intercept' should be a named list.")
    priorEvent_dist_for_intercept <- priorEvent_intercept$dist
    priorEvent_scale_for_intercept <- priorEvent_intercept$scale
    priorEvent_mean_for_intercept <- priorEvent_intercept$location
    priorEvent_df_for_intercept <- priorEvent_intercept$df 
    priorEvent_df_for_intercept[is.na(priorEvent_df_for_intercept)] <- 1
    
    if (!priorEvent_dist_for_intercept %in% unlist(ok_intercept_dists))
      stop("The prior distribution for the event model intercept should be one of ",
           paste(names(ok_intercept_dists), collapse = ", "))
    priorEvent_dist_for_intercept <- 
      ifelse(priorEvent_dist_for_intercept == "normal", 1L, 2L)
    priorEvent_scale_for_intercept <- 
      set_prior_scale(priorEvent_scale_for_intercept, default = 10, 
                      link = family$link)
    priorEvent_df_for_intercept <- min(.Machine$double.xmax, priorEvent_df_for_intercept)
  }

   # Priors for association parameters
  priorAssoc_scaled <- priorAssoc_ops$scaled
  priorAssoc_min_prior_scale <- priorAssoc_ops$min_prior_scale
  
  if (is.null(priorAssoc)) {
    priorAssoc_dist <- 0L
    priorAssoc_mean <- as.array(rep(0, a_K))
    priorAssoc_scale <- priorAssoc_df <- as.array(rep(1, a_K))
  } else {
    if (!is.list(priorAssoc)) 
      stop("'priorAssoc' should be a named list.")
    priorAssoc_dist <- priorAssoc$dist
    priorAssoc_scale <- priorAssoc$scale
    priorAssoc_mean <- priorAssoc$location
    priorAssoc_df <- priorAssoc$df
    priorAssoc_df[is.na(priorAssoc_df)] <- 1
    if (!priorAssoc_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the event model coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "))
    } else if (priorAssoc_dist %in% c("normal", "t")) {
      priorAssoc_dist <- ifelse(priorAssoc_dist == "normal", 1L, 2L)
      priorAssoc_scale <- set_prior_scale(priorAssoc_scale, default = 2.5, 
                                          link = family$link)
    } else {
      priorAssoc_dist <- ifelse(priorAssoc_dist == "hs", 3L, 4L)
    }
    
    priorAssoc_df <- maybe_broadcast(priorAssoc_df, a_K)
    priorAssoc_df <- as.array(pmin(.Machine$double.xmax, priorAssoc_df))
    priorAssoc_mean <- maybe_broadcast(priorAssoc_mean, a_K)
    priorAssoc_mean <- as.array(priorAssoc_mean)
    priorAssoc_scale <- maybe_broadcast(priorAssoc_scale, a_K)
  }
  
  # Minimum scaling of priors for longitudinal submodel(s)
  if (priorLong_scaled && priorLong_dist > 0L) {
    for (j in 1:m) {
      if (y_K[j] > 0L) {
        if (j == 1L) {
          mark_start <- 1 
          mark_end <- y_K[1]
        } else {
          mark_start <- sum(y_K[1:(j-1)]) + 1
          mark_end <- sum(y_K[1:j])
        }
        if (is_gaussian[[j]]) {
          ss <- 2 * sd(y[[j]])
          priorLong_scale[mark_start:mark_end] <- ss * priorLong_scale[mark_start:mark_end]
          priorLong_scale_for_intercept[[j]] <-  ss * priorLong_scale_for_intercept[[j]]
        }
        if (!QR) 
          priorLong_scale[mark_start:mark_end] <- 
            pmax(priorLong_min_prior_scale, priorLong_scale[mark_start:mark_end] / 
                 apply(xtemp[[j]], 2L, FUN = function(x) {
                   num.categories <- length(unique(x))
                   x.scale <- 1
                   if (num.categories == 2) x.scale <- diff(range(x))
                   else if (num.categories > 2) x.scale <- 2 * sd(x)
                   return(x.scale)
                 }))      
      }

    }
  }
  priorLong_scale <- as.array(pmin(.Machine$double.xmax, priorLong_scale))
  priorLong_scale_for_intercept <- 
    min(.Machine$double.xmax, priorLong_scale_for_intercept)

  # Minimum scaling of priors for event submodel
  if (priorEvent_scaled && priorEvent_dist > 0L) {
    priorEvent_scale <- pmax(priorEvent_min_prior_scale, priorEvent_scale / 
                            apply(e_xtemp, 2L, FUN = function(x) {
                              num.categories <- length(unique(x))
                              e.x.scale <- 1
                              if (num.categories == 2) e.x.scale <- diff(range(x))
                              else if (num.categories > 2) e.x.scale <- 2 * sd(x)
                              return(e.x.scale)
                            }))
  }
  priorEvent_scale <- as.array(pmin(.Machine$double.xmax, priorEvent_scale))
  priorEvent_scale_for_intercept <- 
    min(.Machine$double.xmax, priorEvent_scale_for_intercept)

  # Minimum scaling of priors for association parameters    
  if (priorAssoc_dist > 0L) {
    priorAssoc_scale <- pmax(priorAssoc_min_prior_scale, priorAssoc_scale)
  }
  priorAssoc_scale <- as.array(pmin(.Machine$double.xmax, priorAssoc_scale))

  # QR not yet implemented for stan_jm  
  if (QR) {
    stop("QR decomposition is not yet supported by stan_jm or stan_jm.fit")
    if (ncol(xtemp) <= 1)
      stop("'QR' can only be specified when there are multiple predictors.")
    cn <- colnames(xtemp)
    decomposition <- qr(xtemp)
    sqrt_nm1 <- sqrt(nrow(xtemp) - 1L)
    Q <- qr.Q(decomposition)
    R_inv <- qr.solve(decomposition, Q) * sqrt_nm1
    xtemp <- Q * sqrt_nm1
    colnames(xtemp) <- cn
    xbar <- c(xbar %*% R_inv)
  }
 
  #=========================
  # Data for export to Stan
  #=========================

  standata <- list(  
    # dimensions
    m = as.integer(m),
    Npat = as.integer(Npat),
    y_N = as.array(y_N), 
    y_K = as.array(y_K), 
    sum_y_K = as.integer(sum_y_K),
    e_K = as.integer(e_K),
    a_K = as.integer(a_K),
    quadnodes = quadnodes,
    Npat_times_quadnodes = as.integer(Npat * quadnodes),
    sum_y_has_intercept = as.integer(sum_y_has_intercept), 
    
    # data for longitudinal submodel(s)
    family = as.array(family)
    link = as.array(link),
    y_centre = as.array(as.integer(centreLong)),
    y_has_intercept = as.array(y_has_intercept),
    y_has_weights = as.array(y_has_weights),
    y_has_offset = as.array(y_has_offset),   
    y1 = else_empty_vector(y[[1]]),
    y2 = else_empty_vector(y[[2]]),
    y3 = else_empty_vector(y[[3]]),
    y1_xbar = else_empty_vector(xbar[[1]]),
    y2_xbar = else_empty_vector(xbar[[2]]),
    y3_xbar = else_empty_vector(xbar[[3]]),
    y1_X = else_empty_matrix(xtemp[[1]]),
    y2_X = else_empty_matrix(xtemp[[2]]),
    y3_X = else_empty_matrix(xtemp[[3]]),
    y1_weights = else_empty_vector(y_weights[[1]]),
    y2_weights = else_empty_vector(y_weights[[2]]),
    y3_weights = else_empty_vector(y_weights[[3]]),
    y1_offset = else_empty_vector(y_offset[[1]]),
    y2_offset = else_empty_vector(y_offset[[2]]),
    y3_offset = else_empty_vector(y_offset[[3]]),
    
    # data for event submodel
    basehaz_weibull = as.integer(basehaz == "weibull")
    e_centre = as.integer(centreEvent),
    e_has_intercept = as.integer(e_has_intercept),
    nrow_y_Xq = NROW(xqtemp[[1]]),
    nrow_e_Xq = NROW(e_xtemp),
    y1_Xq = else_empty_matrix(xqtemp[[1]]),
    y2_Xq = else_empty_matrix(xqtemp[[2]]),
    y3_Xq = else_empty_matrix(xqtemp[[3]]),
    e_Xq = e_xtemp,
    e_times = c(eventtime, unlist(quadpoint)),
    e_d = c(d, rep(1, length(unlist(quadpoint)))),
    e_xbar = if (centreEvent) as.array(e_xbar) else double(0),
    quadweight_times_half_eventtime = quadweight_times_half_eventtime,
    
    # data for association structure
    assoc = as.integer(assoc),
    has_assoc_ev = as.array(has_assoc_ev),
    has_assoc_es = as.array(has_assoc_es),
    has_assoc_cv = as.array(has_assoc_cv),
    has_assoc_cs = as.array(has_assoc_cs),
    sum_has_assoc_ev = as.integer(sum(has_assoc_ev)),
    sum_has_assoc_es = as.integer(sum(has_assoc_es)),
    sum_has_assoc_cv = as.integer(sum(has_assoc_cv)),
    sum_has_assoc_cs = as.integer(sum(has_assoc_cs)),
    which_b_for_assoc = as.array(which_b_for_assoc),
    size_which_b_for_assoc = as.integer(length(size_which_b_for_assoc)),
    
    # priors
    priorLong_dist = priorLong_dist, 
    priorLong_dist_for_intercept = priorLong_dist_for_intercept,  
    priorEvent_dist = priorEvent_dist,
    priorEvent_dist_for_intercept = priorEvent_dist_for_intercept,
    priorAssoc_dist = priorAssoc_dist,    
    
    # hyperparameters for priors
    priorLong_mean = priorLong_mean, 
    priorLong_mean_for_intercept = priorLong_mean_for_intercept,
    priorEvent_mean = priorEvent_mean, 
    priorEvent_mean_for_intercept = priorEvent_mean_for_intercept,
    priorAssoc_mean = priorAssoc_mean, 
    priorLong_scale = priorLong_scale, 
    priorLong_scale_for_intercept = priorLong_scale_for_intercept, 
    priorEvent_scale = priorEvent_scale, 
    priorEvent_scale_for_intercept = priorEvent_scale_for_intercept, 
    priorAssoc_scale = priorAssoc_scale, 
    priorLong_df = priorLong_df, 
    priorLong_df_for_intercept = priorLong_df_for_intercept,  
    priorEvent_df = priorEvent_df, 
    priorEvent_df_for_intercept = priorEvent_df_for_intercept,
    priorAssoc_df = priorAssoc_df, 
    priorLong_scale_for_dispersion = priorLong_scale_for_dispersion,
    priorEvent_scale_for_weibull = priorEvent_scale_for_weibull,
    
    prior_PD = as.integer(prior_PD)
  )  
  
  # data for random effects
  group <- lapply(seq_along(m), function(x) {
                    pad_reTrms(Z = Z[[x]], 
                               cnms = cnms[[x]], 
                               flist = flist[[x]])})
  Z <- lapply(seq_along(group), function(x) group[[x]]$Z)
  Z <- lapply(Z, else_empty_matrix)
  cnms <- lapply(seq_along(group), function(x) group[[x]]$cnms)
  flist <- lapply(seq_along(group), function(x) group[[x]]$flist)
  p_y <- sapply(cnms, function(x) sapply(x, length))  # ranefs in each submodel
  p <- sum(p_y)  # total number of ranefs
  l <- sapply(attr(flist[[1]], "assign"), function(i) 
    nlevels(group$flist[[i]]))
  #t <- length(p)  # set to 1 for stan_jm 
  #######! Needs changing since b is not ordered by submodel, but ordered by patient
  group_nms <- names(cnms[[1]])
  b_nms <- character()
  for (x in 1:m) {
    for (i in seq_along(cnms[[x]])) {
      # if you change this change .pp_data_mer_z() as well
      nm <- group_nms[i]
      nms_i <- paste(cnms[[x]][[i]], nm)
      if (length(nms_i) == 1) {
        b_nms <- c(b_nms, paste0(nms_i, ":", levels(flist[[x]][[nm]])))
      } else {
        b_nms <- c(b_nms, c(t(sapply(nms_i, paste0, ":", levels(flist[[nm]])))))
      }
    }
  }
  ########  
  g_nms <- unlist(lapply(1:m, FUN = function(i) {
    paste0(cnms[[i]][[1]], paste0("|Submodel ", i))
  }))
  standata$t <- 1
  standata$p_y <- as.array(p_y)    
  standata$p <- as.array(p)
  standata$l <- as.array(l)
  standata$q <- sum(sapply(Z, ncol))
  standata$q_y1 <- ncol(Z[[1]])
  standata$q_y2 <- ncol(Z[[2]])
  standata$q_y3 <- ncol(Z[[3]])
  standata$len_theta_L <- sum(choose(p, 2), p)
  parts <- lapply(Z, extract_sparse_parts)
  standata$num_non_zero <- as.array(sapply(1:3, function(x) length(parts[[x]]$w)))
  standata$w1 <- parts[[1]]$w
  standata$w2 <- parts[[2]]$w
  standata$w3 <- parts[[3]]$w
  standata$v1 <- parts[[1]]$v
  standata$v2 <- parts[[2]]$v
  standata$v3 <- parts[[3]]$v
  standata$u1 <- parts[[1]]$u
  standata$u2 <- parts[[2]]$u
  standata$u3 <- parts[[3]]$u
 
  # data for random effects in GK quadrature
  Z <- lapply(Z, else_empty_matrix)
  parts_Zq <- lapply(Zq, extract_sparse_parts)
  standata$num_non_zero_Zq <- as.array(sapply(1:3, function(x) length(parts_Zq[[x]]$w)))
  standata$w1_Zq <- parts_Zq[[1]]$w
  standata$w2_Zq <- parts_Zq[[2]]$w
  standata$w3_Zq <- parts_Zq[[3]]$w
  standata$v1_Zq <- parts_Zq[[1]]$v
  standata$v2_Zq <- parts_Zq[[2]]$v
  standata$v3_Zq <- parts_Zq[[3]]$v
  standata$u1_Zq <- parts_Zq[[1]]$u
  standata$u2_Zq <- parts_Zq[[2]]$u
  standata$u3_Zq <- parts_Zq[[3]]$u

  # hyperparameters for random effects model
  standata$shape <- as.array(maybe_broadcast(decov$shape, t))
  standata$scale <- as.array(maybe_broadcast(decov$scale, t))
  standata$len_concentration <- sum(p[p > 1])
  standata$concentration <- 
    as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
  standata$len_regularization <- sum(p > 1)
  standata$regularization <- 
    as.array(maybe_broadcast(decov$regularization, sum(p > 1))) 
    
  standata$family <- switch(family$family, 
                            gaussian = 1L, 
                            Gamma = 2L,
                            inverse.gaussian = 3L,
                            binomial = 4L,
                            poisson = 5L,
                            "neg_binomial_2" = 6L)     
  
  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$jm
  pars <- c(if (y_has_intercept[[1]]) "y1_gamma",
            if (y_has_intercept[[2]]) "y2_gamma",
            if (y_has_intercept[[3]]) "y3_gamma",
            "alpha", 
            "y_beta",
            if (y_K[1] > 0) "y1_beta",
            if (y_K[2] > 0) "y2_beta",
            if (y_K[3] > 0) "y3_beta",
            "e_beta",
            if (assoc) "a_beta",
            if (length(group)) "b",
            "y_dispersion", 
            if (standata$basehaz_weibull) "weibull_shape",
            "mean_PPD")
  if (algorithm == "optimizing") {
    out <- optimizing(stanfit, data = standata, 
                      draws = 1000, constrained = TRUE, ...)
    new_names <- names(out$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    if (QR) {
      out$par[mark] <- R_inv %*% out$par[mark]
      out$theta_tilde[,mark] <- out$theta_tilde[, mark] %*% t(R_inv)
    }
    new_names[mark] <- colnames(xtemp)
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    new_names[grepl("dispersion(\\[1\\])?$", new_names)] <- 
      if (is_gaussian) "sigma" else
        if (is_gamma) "shape" else
          if (is_ig) "lambda" else 
            if (is_nb) "overdispersion" else NA
    names(out$par) <- new_names
    colnames(out$theta_tilde) <- new_names
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, 
                                             chains = 0))
    return(out)
    
  } else {
    if (algorithm == "sampling") {
      sampling_args <- set_sampling_args(
        object = stanfit, 
        prior = prior, 
        user_dots = list(...), 
        user_adapt_delta = adapt_delta, 
        data = standata, 
        pars = pars, 
        show_messages = FALSE)
      stanfit <- do.call(sampling, sampling_args)
    } else {
      # meanfield or fullrank vb
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, ...)
      if (algorithm == "meanfield" && !QR) 
        msg_meanfieldQR()
    }
    if (QR) {
      thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, 
                        permuted = FALSE)
      betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
      end <- tail(dim(betas), 1L)
      for (chain in 1:end) for (param in 1:nrow(betas)) {
        stanfit@sim$samples[[chain]][[has_intercept + param]] <-
          if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
      }
    }
    new_names <- c(if (has_intercept) "(Intercept)", 
                   colnames(xtemp), 
                   if (length(group)) c(paste0("b[", b_nms, "]")),
                   if (is_gaussian) "sigma", 
                   if (is_gamma) "shape", 
                   if (is_ig) "lambda",
                   if (is_nb) "overdispersion", 
                   "mean_PPD", 
                   "log-posterior")
    stanfit@sim$fnames_oi <- new_names
    return(stanfit)
  }
}



  Z <- pad_reTrms(Z = t(group$Zt), cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call, terms = NULL, model = NULL, 
               prior.info = get_prior_info(call, formals()),
               na.action, contrasts, algorithm, mod)
  out <- stanreg(fit)
  
  return(out)
}


# If object is NA or NULL then return empty vector or matrix
#
# @param x Object to check
else_empty_vector <- function(x) {
  if (!is.na(x) && !is.null(x)) {
    x
  } else {
    double()
  }
}
else_empty_matrix <- function(x) {
  if (!is.na(x) && !is.null(x)) {
    x
  } else {
    matrix(,0,0)
  } 
}

# Extend vector or list to length n replacing the
# additional new elements with integer or logical 'new'
#
# @param x Vector or list to extend
# @param n Integer specifying the desired length of
#   the returned vector or list
# @param new The entry that the additional elements
#   should be set to
maybe_extend <- function(x, n, new) {
  if (!length(x)) {
    rep(new, times = n)
  } else if (length(x) < n) {
    c(x, rep(new, times = (n - length(x))))
  } else {
    x
  }
}

# Add extra level _NEW_ to each group
# 
# @param Z ranef indicator matrix
# @param cnms group$cnms
# @param flist group$flist
pad_reTrms <- function(Z, cnms, flist) {
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  last <- cumsum(l * p)
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])), 
                            paste0("_NEW_", names(flist)[i]))
  }
  n <- nrow(Z)
  mark <- length(p) - 1L
  if (getRversion() < "3.2.0") {
    Z <- cBind(Z, Matrix(0, nrow = n, ncol = p[length(p)], sparse = TRUE))
    for (i in rev(head(last, -1))) {
      Z <- cBind(cBind(Z[, 1:i, drop = FALSE],
                       Matrix(0, n, p[mark], sparse = TRUE)),
                 Z[, (i+1):ncol(Z), drop = FALSE])
      mark <- mark - 1L
    }
  }
  else {
    Z <- cbind2(Z, Matrix(0, nrow = n, ncol = p[length(p)], sparse = TRUE))
    for (i in rev(head(last, -1))) {
      Z <- cbind(Z[, 1:i, drop = FALSE],
                 Matrix(0, n, p[mark], sparse = TRUE),
                 Z[, (i+1):ncol(Z), drop = FALSE])
      mark <- mark - 1L
    }
  }
  nlist(Z, cnms, flist)
}

# Drop the extra reTrms from a matrix x
#
# @param x A matrix (e.g. the posterior sample or matrix of summary stats)
# @param columns Do the columns (TRUE) or rows (FALSE) correspond to the
#   variables?
unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x))
    return(unpad_reTrms.matrix(x, ...))
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}
unpad_reTrms.matrix <- function(x, columns = TRUE, ...) {
  nms <- if (columns) 
    colnames(x) else rownames(x)
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (columns) x[, keep, drop = FALSE] else x[keep, , drop = FALSE]
}
 