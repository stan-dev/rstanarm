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

#' Bayesian joint longitudinal and time-to-event models via Stan
#' 
#' Fits a shared parameter joint model for longitudinal and time-to-event 
#' (e.g. survival) data under a Bayesian framework using Stan.
#' 
#' @param formulaLong A two-sided linear formula object describing both the 
#'   fixed-effects and random-effects parts of the longitudinal submodel  
#'   (see \code{\link[lme4]{glmer}} for details). For a multivariate joint 
#'   model (i.e. more than one longitudinal marker) this should 
#'   be a list of such formula objects, with each element
#'   of the list providing the formula for one of the longitudinal submodels.
#' @param dataLong A data frame containing the variables specified in
#'   \code{formulaLong}. If fitting a multivariate joint model, then this can
#'   be either a single data frame which contains the data/variables for all 
#'   the longitudinal submodels, or it can be a list of data frames where each
#'   element of the list provides the data for one of the longitudinal 
#'   submodels.
#' @param formulaEvent A two-sided formula object describing the event
#'   submodel. The left hand side of the formula should be a \code{Surv()} 
#'   object. See \code{\link[survival]{Surv}}.
#' @param dataEvent A data frame containing the variables specified in
#'   \code{formulaEvent}.
#' @param time_var A character string specifying the name of the variable 
#'   in \code{dataLong} which represents time.
#' @param id_var A character string specifying the name of the variable in
#'   \code{dataLong} which distinguishes between individuals. This can be
#'   left unspecified if there is only one grouping factor (which is assumed
#'   to be the individual). If there is more than one grouping factor (i.e.
#'   clustering beyond the level of the individual) then the \code{id_var}
#'   argument must be specified.
#' @param family The family (and possibly also the link function) for the 
#'   longitudinal submodel(s). See \code{\link[lme4]{glmer}} for details. 
#'   If fitting a multivariate joint model, then this can optionally be a
#'   list of families, in which case each element of the list specifies the
#'   family for one of the longitudinal submodels.
#' @param assoc A character string or character vector specifying the joint
#'   model association structure. Possible association structures that can
#'   be used include: "etavalue" (the default); "etaslope"; "etalag"; "etaauc; "muvalue"; 
#'   "muslope"; "mulag"; "muauc"; "shared_b"; "shared_coef"; or "null". 
#'   These are described in the 
#'   \strong{Details} section below. For a multivariate joint model, 
#'   different association structures can optionally be used for 
#'   each longitudinal submodel by specifying a list of character
#'   vectors, with each element of the list specifying the desired association 
#'   structure for one of the longitudinal submodels. Specifying \code{assoc = NULL}
#'   will fit a joint model with no association structure (equivalent  
#'   to fitting separate longitudinal and time-to-event models). Also see the
#'   \strong{Examples} section below.
#' @param base_haz A character string indicating which baseline hazard to use
#'   for the event submodel. Options are a Weibull baseline hazard
#'   (\code{"weibull"}, the default), a B-splines approximation estimated 
#'   for the log baseline hazard (\code{"bs"}), or a piecewise
#'   constant baseline hazard (\code{"piecewise"}).
#' @param df An optional positive integer specifying the degrees of freedom 
#'   for the B-splines if \code{base_haz = "bs"}, or the number of
#'   intervals used for the piecewise constant baseline hazard if 
#'   \code{base_haz = "piecewise"}. The default is 6.
#' @param knots An optional numeric vector specifying the internal knot 
#'   locations for the B-splines if \code{base_haz = "bs"}, or the 
#'   internal cut-points for defining intervals of the piecewise constant 
#'   baseline hazard if \code{base_haz = "piecewise"}. Knots cannot be
#'   specified if \code{df} is specified. If not specified, then the 
#'   default is to use \code{df - 4} knots if \code{base_haz = "bs"},
#'   or \code{df - 1} knots if \code{base_haz = "piecewise"}, which are
#'   placed at equally spaced percentiles of the distribution of
#'   observed event times.
#' @param quadnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard in the likelihood function. 
#'   Options are 15 (the default), 11 or 7.
#' @param subsetLong,subsetEvent Same as subset in \code{\link[stats]{glm}}.
#'   However, if fitting a multivariate joint model and a list of data frames 
#'   is provided in \code{dataLong} then a corresponding list of subsets 
#'   must be provided in \code{subsetLong}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param weights Experimental and should be used with caution. 
#'   The user can optionally supply a 2-column data frame containing a set of
#'   'prior weights' to be used in the estimation process. The data frame should
#'   contain two columns: the first containing the IDs for each individual, and 
#'   the second containing the corresponding weights. The data frame should only
#'   have one row for each individual; that is, weights should be constant 
#'   within individuals.
#' @param offset Not currently implemented. Same as \code{\link[stats]{glm}}. 
#' @param centreLong,centreEvent A logical specifying whether the predictor
#'   matrix for the longitudinal submodel(s) or event submodel should be 
#'   centred. 
#' @param init The method for generating the initial values for the MCMC.
#'   The default is \code{"model_based"}, which uses those obtained from 
#'   fitting separate longitudinal and time-to-event models prior  
#'   to fitting the joint model. Parameters that cannot be obtained from 
#'   fitting separate longitudinal and time-to-event models are initialised 
#'   at 0. This provides reasonable initial values which should aid the MCMC
#'   sampler. However, it is recommended that any final analysis should be
#'   performed with several MCMC chains each initiated from a different
#'   set of initial values; this can be obtained by setting
#'   \code{init = "random"}. Other possibilities for specifying \code{init}
#'   are those described for \code{\link[rstan]{stan}}.     
#' @param ... Further arguments passed to 
#'   \code{\link[rstan]{sampling}} (e.g. \code{iter}, \code{chains}, 
#'   \code{cores}, etc.).
#' @param priorLong,priorEvent,priorAssoc The prior distributions for the 
#'   regression coefficients in the longitudinal submodel(s), event submodel,
#'   and the association parameter(s). 
#'   Can be a call to \code{normal}, \code{student_t},
#'   \code{cauchy}, \code{hs} or \code{hs_plus}. See \code{\link{priors}} for
#'   details. To to omit a prior --- i.e., to use a flat (improper) uniform
#'   prior --- set equal to \code{NULL}.
#' @param priorLong_intercept,priorEvent_intercept The prior distributions  
#'   for the intercepts in the longitudinal submodel(s) and event submodel.
#'   Can be a call to \code{normal}, \code{student_t} or
#'   \code{cauchy}. See \code{\link{priors}} for details. To to omit a prior
#'   --- i.e., to use a flat (improper) uniform prior --- set
#'   equal to \code{NULL}. (\strong{Note:} If \code{centreLong} or 
#'   \code{centreEvent} is set to \code{TRUE} then the prior
#'   distribution for the intercept is set so it applies to the value when all
#'   predictors are centered).
#' @param priorLong_ops,priorEvent_ops,priorAssoc_ops Additional options 
#'   related to prior distributions for the longitudinal submodel(s),
#'   event submodel, and the association parameter(s). To omit a prior on the 
#'   dispersion parameters in the longitudinal submodel(s) set 
#'   \code{priorLong_ops} to \code{NULL}. To omit a prior on the 
#'   Weibull shape parameter (if \code{base_haz = "weibull"}) or the
#'   spline coefficients (if \code{base_haz = "bs"}) set 
#'   \code{priorEvent_ops} to \code{NULL}. Otherwise see \code{\link{priors}}. 
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#' @param prior_PD A logical (defaulting to \code{FALSE}) indicating
#'   whether to draw from the prior predictive distribution instead of
#'   conditioning on the data.
#' @param adapt_delta See \code{\link{adapt_delta}} for details.
#' @param max_treedepth A positive integer specifying the maximum treedepth 
#'   for the non-U-turn sampler. See the \code{control} argument in 
#'   \code{\link[rstan]{stan}}.
#' @param QR QR decomposition of the predictor matrix. Not yet implemented.
#' @param algorithm Character string (possibly abbreviated) indicating the 
#'   estimation approach to use. Currently, only "sampling" (for MCMC) is 
#'   allowed.
#' @param debug Should not be used. This is used in package development to
#'   help identify the location of bugs in the Stan code. When set to 
#'   \code{TRUE}, the Stan file executes several print statements which 
#'   can help to identify the location of bugs or numerical instabilities.  
#'   
#' @details The \code{stan_jm} function can be used to fit a joint model (also 
#'   known as a shared parameter model) for longitudinal and time-to-event data 
#'   under a Bayesian framework. 
#'   The joint model may be univariate (with only one longitudinal submodel) or
#'   multivariate (with more than one longitudinal submodel). Multi-level 
#'   clustered data are allowed (e.g. patients within clinics), provided that the
#'   individual (e.g. patient) is the lowest level of clustering. The underlying
#'   estimation is carried out using the Bayesian C++ package Stan 
#'   (\url{http://mc-stan.org/}). \cr
#'   \cr 
#'   For the longitudinal submodel a generalised linear mixed model is assumed 
#'   with any of the \code{\link[stats]{family}} choices allowed by 
#'   \code{\link[lme4]{glmer}}. \cr
#'   \cr
#'   For the event submodel a parametric
#'   proportional hazards model is assumed. The baseline hazard can be estimated 
#'   using either a Weibull distribution (\code{base_haz = "weibull"}) or a
#'   piecewise constant baseline hazard (\code{base_haz = "piecewise"}), or 
#'   approximated using cubic B-splines (\code{base_haz = "bs"}). 
#'   If either of the latter two are used then the degrees of freedom, 
#'   or the internal knot locations, can be optionally specified. If
#'   the degrees of freedom are specified (through the \code{df} argument) then
#'   the knot locations are automatically generated based on the 
#'   distribution of the observed event times (not including censoring times). 
#'   Otherwise internal knot locations can be specified 
#'   directly through the \code{knots} argument. If neither \code{df} or
#'   \code{knots} is specified, then the default is to set \code{df} equal to 6.
#'   It is not possible to specify both \code{df} and \code{knots}. \cr
#'   \cr
#'   The association structure for the joint model can be based on any of the 
#'   following parameterisations: current value of the linear predictor in the 
#'   longitudinal submodel (\code{"etavalue"}); first derivative 
#'   (slope) of the linear predictor in the longitudinal submodel 
#'   (\code{"etaslope"}); lagged value of the linear predictor in the 
#'   longitudinal submodel (\code{"etalag(#)"}, replacing \code{#} with the
#'   desired lag in units of the time variable); the area under the curve of 
#'   the linear predictor in the longitudinal submodel (\code{"etaauc"}); 
#'   current expected value of the 
#'   longitudinal submodel (\code{"muvalue"}); lagged expected value of the
#'   longitudinal submodel (\code{"mulag(#)"}, replacing \code{#} with the
#'   desired lag in units of the time variable); the area under the curve of 
#'   the expected value from the longitudinal submodel (\code{"muauc"}); 
#'   shared individual-level random  
#'   effects (\code{"shared_b"}); shared individual-level random effects which also
#'   incorporate the corresponding fixed effect as well as any corresponding 
#'   random effects for clustering levels higher than the individual)
#'   (\code{"shared_coef"}); or no 
#'   association structure (equivalent to fitting separate longitudinal 
#'   and event models) (\code{"null"} or \code{NULL}). 
#'   More than one association structure can be specified, however,
#'   not all possible combinations are allowed.   
#'   Note that for the lagged association structures (\code{"etalag(#)"} and 
#'   \code{"mulag(#)"}) use baseline values (time = 0) for the instances where the 
#'   time lag results in a time prior to baseline. When using the 
#'   \code{"etaauc"} or \code{"muauc"} association structures, the area under
#'   the curve is evaluated using Gauss-Kronrod quadrature with 15 quadrature nodes. 
#'   By default, \code{"shared_b"} and \code{"shared_coef"} contribute all 
#'   random effects to the association structure; however, a subset of the random effects can 
#'   be chosen by specifying their indices between parentheses as a suffix, for 
#'   example, "shared_b(1)" or "shared_b(1:3)" or "shared_b(1,2,4)", and so on. \cr
#'   \cr
#'   Time-varying covariates are allowed in both the 
#'   longitudinal and event submodels. These should be specified in the data 
#'   in the same way as they normally would when fitting a separate 
#'   longitudinal model using \code{\link[lme4]{lmer}} or a separate 
#'   time-to-event model using \code{\link[survival]{coxph}}. These time-varying
#'   covariates should be exogenous in nature, otherwise they would perhaps 
#'   be better specified as an additional outcome (i.e. by including them as an 
#'   additional longitudinal outcome/submodel in the joint model). \cr
#'   \cr
#'   Bayesian estimation of the joint model is performed via MCMC. The Bayesian  
#'   model includes independent priors on the 
#'   regression coefficients for both the longitudinal and event submodels, 
#'   including the association parameter(s) (in much the same way as the
#'   regression parameters in \code{\link{stan_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters (in the same way as \code{\link{stan_glmer}}). 
#'   See \code{\link{priors}} for more information about the priors distributions
#'   that are available. \cr
#'   \cr
#'   Gauss-Kronrod quadrature is used to numerically evaluate the integral  
#'   over the cumulative hazard in the likelihood function for the event submodel.
#'   The accuracy of the numerical approximation can be controlled using the
#'   number of quadrature nodes, specified through the \code{quadnodes} 
#'   argument. Using a higher number of quadrature nodes will result in a more 
#'   accurate approximation.
#'
#' @return A \link[=stanjm-object]{stanjm} object is returned.
#' 
#' @template reference-stan-manual
#' 
#' @seealso \code{\link{stanjm-object}}, \code{\link{stanjm-methods}}, 
#'   \code{\link{print.stanjm}}, \code{\link{summary.stanjm}},
#'   \code{\link{posterior_traj}}, \code{\link{posterior_survfit}}, 
#'   \code{\link{posterior_predict}}, \code{\link{posterior_interval}},
#'   \code{\link{pp_check}}, \code{\link{ps_check}}.
#'    
#' @examples
#' \donttest{
#' #####
#' # Univariate joint model, with association structure based on the 
#' # current value of the linear predictor
#' f1 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong_subset,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv_subset,
#'               time_var = "year")
#' summary(f1) 
#'         
#' #####
#' # Univariate joint model, with association structure based on the 
#' # current value of the linear predictor and shared random intercept
#' f2 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong_subset,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv_subset,
#'               assoc = c("etavalue", "shared_b"),
#'               time_var = "year")
#' summary(f2)          
#' 
#' ######
#' # Multivariate joint model, with association structure based 
#' # on the current value of the linear predictor in each longitudinal 
#' # submodel and shared random intercept from the second longitudinal 
#' # submodel only (which is the first random effect in that submodel
#' # and is therefore indexed the '(1)' suffix in the code below)
#' mv1 <- stan_jm(
#'         formulaLong = list(
#'           logBili ~ year + (1 | id), 
#'           albumin ~ sex + year + (1 + year | id)),
#'         dataLong = pbcLong_subset,
#'         formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'         dataEvent = pbcSurv_subset,
#'         assoc = list("etavalue", c("etavalue", "shared_b(1)")), 
#'         time_var = "year")
#' summary(mv1)
#' 
#' # To include both the random intercept and random slope in the shared 
#' # random effects association structure for the second longitudinal 
#' # submodel, we could specify the following:
#' #   update(mv1, assoc = list("etavalue", c("etavalue", "shared_b"))
#' # which would be equivalent to:  
#' #   update(mv1, assoc = list("etavalue", c("etavalue", "shared_b(1,2)"))
#' # or:
#' #   update(mv1, assoc = list("etavalue", c("etavalue", "shared_b(1:2)"))     
#' 
#' ######
#' # Multivariate joint model, estimated using multiple MCMC chains 
#' # run in parallel across all available PC cores
#' mv2 <- stan_jm(formulaLong = list(
#'         logBili ~ year + (1 | id), 
#'         albumin ~ sex + year + (1 +  year | id)),
#'         dataLong = pbcLong_subset,
#'         formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'         dataEvent = pbcSurv_subset,
#'         assoc = list("etavalue", c("etavalue", "shared_b(1)")),
#'         time_var = "year",
#'         chains = 3, refresh = 25,
#'         cores = parallel::detectCores())
#' summary(mv2)  
#' 
#' #####
#' # Here we provide an example of specifying an association structure 
#' # based on the lagged value of the linear predictor, where the lag
#' # is 2 time units (i.e. 2 years in this example)
#' f3 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong_subset,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv_subset,
#'               time_var = "year",
#'               assoc = "etalag(2)")
#' summary(f3) 
#' 
#' #####
#' # Here we provide an example of specifying an association structure with 
#' # interaction terms. Here we specify that we want to use an association
#' # structure based on the current value of the linear predictor from
#' # the longitudinal submodel ("etavalue"), but we will also specify
#' # that we want to interact this with the treatment covariate (trt) from
#' # pbcLong_subset data frame so that we can estimate a different association 
#' # parameter (i.e. estimated effect of log serum bilirubin on the log hazard 
#' # of death) for each treatment group
#' f4 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong_subset,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv_subset,
#'               time_var = "year",
#'               assoc = c("etavalue", "etavalue_interact(~ trt)"))
#' summary(f4)  
#'                     
#' }
#' 
#' @export
#' @import data.table
#' @importFrom survival Surv
#' 
stan_jm <- function(formulaLong, dataLong, 
                    formulaEvent, dataEvent, 
                    time_var, id_var, family = gaussian,
                    assoc = "etavalue", interaction_data,
                    base_haz = c("weibull", "bs", "piecewise"), 
                    df, knots, quadnodes = 15, 
                    subsetLong, subsetEvent, 
                    na.action = getOption("na.action", "na.omit"),
                    weights, offset, contrasts,
                    centreLong = FALSE, centreEvent = FALSE, 
                    init = "model_based", ...,				          
                    priorLong = normal(), priorLong_intercept = normal(),
                    priorLong_ops = priorLong_options(),
                    priorEvent = normal(), priorEvent_intercept = normal(),
                    priorEvent_ops = priorEvent_options(),
                    priorAssoc = normal(),
                    priorAssoc_ops = priorAssoc_options(),
                    prior_covariance = decov(), prior_PD = FALSE,
                    adapt_delta = NULL, max_treedepth = NULL, QR = FALSE, 
                    algorithm = c("sampling", "meanfield", "fullrank")) {
  
  
  #=============================
  # Pre-processing of arguments
  #=============================  
  
  # Check for arguments not yet implemented
  if (!missing(offset)) 
    stop("Offsets are not yet implemented for stan_jm")
  if (QR) 
    stop("QR decomposition not yet implemented for stan_jm")   
  algorithm <- match.arg(algorithm)
#  if (algorithm %in% c("meanfield", "fullrank"))
#    stop ("Meanfield and fullrank algorithms not yet implemented for stan_jm")
#  if ((init == "model_based") && any(unlist(c(centreLong, centreEvent)))) 
#    stop("Cannot use model based initial values when 'centreLong = TRUE'",
#         " or 'centreEvent = TRUE'.")  
  if (missing(weights)) weights <- NULL
  if (missing(id_var)) id_var <- NULL
  if (missing(subsetLong)) subsetLong <- NULL
  if (missing(interaction_data)) interaction_data <- NULL
  
  # Matched call
  call <- match.call(expand.dots = TRUE)    
  mc <- match.call(expand.dots = FALSE)
  mc$time_var <- mc$id_var <- mc$assoc <- mc$base_haz <- 
    mc$df <- mc$knots <- mc$quadnodes <- 
    mc$centreLong <- mc$centreEvent <- NULL
  mc$priorLong <- mc$priorLong_intercept <- mc$priorLong_ops <- 
    mc$priorEvent <- mc$priorEvent_intercept <- mc$priorEvent_ops <-
    mc$priorAssoc <- mc$priorAssoc_ops <- mc$prior_covariance <-
    mc$prior_PD <- mc$algorithm <- mc$scale <- 
    mc$concentration <- mc$shape <- mc$init <- 
    mc$adapt_delta <- mc$max_treedepth <- 
    mc$... <- mc$QR <- NULL
  mc$weights <- NULL
  
  # Create call for longitudinal submodel  
  y_mc <- mc
  y_mc <- strip_nms(y_mc, "Long") 
  y_mc$formulaEvent <- y_mc$dataEvent <- y_mc$subsetEvent <- NULL
  
  # Is formulaLong a list?
  if (is(eval(y_mc$formula), "formula")) { # number of long. markers
    formula_list <- FALSE
    M <- 1L
  } else if (is.list(eval(y_mc$formula))) {
    formula_list <- TRUE
    M <- length(eval(y_mc$formula))
    if (!all(sapply(eval(y_mc$formula), function(x)
      (is(x, "formula")))))
      stop("'formulaLong' should be a formula object or a list of formula ",
           "objects")               
  } else {
    stop("'formulaLong' should be a formula object or, for a multivariate ",
         "joint model, a list of formula objects with length equal to the ",
         "desired number of longitudinal markers")
  }
  
  # Is dataLong a list?
  if (is.data.frame(eval(y_mc$data))) {
    data_list <- FALSE
  } else if (is(eval(y_mc$data), "list")) {
    data_list <- TRUE
    if (length(eval(y_mc$data)) != M)
      stop("dataLong appears to be a list of the incorrect length")
  } else {
    stop("'dataLong' should be a data frame or possibly a list of data ",
         "frames. The latter is only required when fitting a multivariate ",
         "joint model using different data for each longitudinal submodel.")
  }
  
  # Is subset a list?
  if (is.null(y_mc$subset)) {
    y_subset_list <- NULL
  } else if (is(eval(y_mc$subset), "list")) {
    y_subset_list <- TRUE
    if (length(eval(y_mc$subset)) != M)
      stop("subsetLong appears to be a list of the incorrect length")
  } else if (is.vector(eval(y_mc$subset))) {
    y_subset_list <- FALSE
  } else {
    stop("'subsetLong' should be a vector or possibly a list of vectors. ",
         "The latter is only required if fitting a multivariate joint ",
         "model and using a different subset of data for each ",
         "longitudinal submodel.")
  }
  
  # Is family a list?
  if (is.null(eval(y_mc$family))) {
    family_list <- NULL
  } else if (is(eval(y_mc$family), "list")) {
    family_list <- TRUE
    if (length(eval(y_mc$family)) != M)
      stop("'family' argument must be a family or possibly a list of families. ",
           "The latter is only required when fitting a multivariate joint ",
           "model with a different family and/or link function for some ",
           "of the longitudinal submodels.")
  } else family_list <- FALSE
  if (is.null(family_list) || !family_list) family <- list(family)  # convert to list
  if (length(family) < M) family <- rep(family, M)  # repeat family if necessary
  
  # Is assoc a list? If not, then convert to list
  if (is.list(assoc)) {  # if list, then check length
    if (!(length(assoc) %in% c(1,M)))
      stop("`assoc' should be a list of length 1 or length equal to the ",
           "number of longitudinal markers")
  } else if (is.character(assoc) || is.null(assoc)) {  # if not list, then convert to list
    assoc <- list(assoc)
  } else {  # else return error
    stop("'assoc' should be NULL, a character string, a character vector or, for a ",
         "multivariate joint model, possibly a list of character strings ",
         "or character vectors. The latter is only required if using a different ",
         "association structure for linking each longitudinal submodel to the ",
         "event outcome.")
  }
  if (length(assoc) != M) assoc <- rep(assoc, M)
  
  # Check family and link
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  family <- lapply(family, validate_family)
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
                 family, seq_along(family), SIMPLIFY = TRUE)
  if (any(lapply(link, length) == 0L)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  # Check for weights
  has_weights <- (!is.null(weights))
  
  # Create call for each longitudinal submodel separately
  m_mc <- list()  # list containing matched calls for each marker
  for (m in 1:M) {
    m_mc[[m]] <- y_mc
    m_mc[[m]]$formula <- if (formula_list) eval(y_mc$formula)[[m]] else y_mc$formula
    m_mc[[m]]$data    <- if (data_list)    eval(y_mc$data)[[m]]    else y_mc$data
    if (!is.null(y_subset_list))   
      m_mc[[m]]$subset  <- if (y_subset_list) eval(y_mc$subset)[[m]]  else y_mc$subset
    if (!is.null(family_list))   
      m_mc[[m]]$family  <- if (family_list)   eval(y_mc$family)[[m]]  else y_mc$family    
  }
  
  # Create call for event submodel
  e_mc <- mc
  e_mc <- strip_nms(e_mc, "Event")
  e_mc$formulaLong <- e_mc$dataLong <- e_mc$family <-
    e_mc$subsetLong <- NULL
  
  # Standardised GK quadrature points and weights
  quadpoints <- get_quadpoints(quadnodes)
  
  
  #================================
  # Data for longitudinal submodel
  #================================
  
  # Items to store for each longitudinal submodel
  y_mod         <- list()     # fitted long. submodels
  y_vars        <- list()     # the variables used in fitting the model
  y_is_real     <- c()        # indicator of response vector being reals (not integers)
  y             <- list()     # response vector
  x             <- list()     # design matrix with intercept
  xtemp         <- list()     # design matrix without intercept, possibly centred
  xbar          <- list()     # means of predictors
  trials        <- list()     # num. of trials for binomial outcome
  y_centre      <- c()        # submodel has intercept
  y_has_intercept <- c()      # submodel has intercept
  y_has_intercept_unbound <- c()    # has unbounded intercept
  y_has_intercept_lobound <- c()    # has lower bounded intercept
  y_has_intercept_upbound <- c()    # has upper bounded intercept
  y_has_dispersion <- c()     # submodel has dispersion term
  y_N           <- c()        # num. observations
  y_N01         <- list()     # num. 0 and 1 observations if bernoulli
  y_real_N      <- c()        # num. observations, for real outcomes
  y_int_N       <- c()        # num. observations, for integer outcomes
  y_K           <- c()        # num. predictors (excluding intercept)
  y_weights     <- list()     # prior weights
  y_offset      <- list()     # offsets
  Ztlist        <- list()     # Z matrices
  y_cnms          <- list()   
  y_flist         <- list()   
  y_gamma         <- c()      # initial values for intercepts
  se_y_gamma      <- c()      
  y_gamma_unbound <- list()   # initial values for intercepts
  y_gamma_lobound <- list()   # initial values for intercepts
  y_gamma_upbound <- list()   # initial values for intercepts
  y_beta          <- list()   # initial values for coefs
  se_y_beta       <- list()  
  y_dispersion    <- c()      # initial values for dispersion
  sd_b            <- list()   # initial values for random effect sds
  b_Cov           <- list()   # initial values for correlation matrix
  b_Corr          <- list()   # initial values for correlation matrix
  ord             <- list()   # ordering of y vector, only present if outcome is bernoulli
  
  for (m in 1:M) {
    
    # Fit separate longitudinal model
    if ((family[[m]]$family == "gaussian") && (family[[m]]$link == "identity")) {
      m_mc[[m]][[1]] <- quote(lme4::lmer)
      m_mc[[m]]$family <- NULL
      m_mc[[m]]$control <- get_control_args()
    } else if (family[[m]]$family == "neg_binomial_2") {
      m_mc[[m]][[1]] <- quote(lme4::glmer.nb)
      m_mc[[m]]$family <- NULL
      m_mc[[m]]$control <- get_control_args(glmer = TRUE)               
    } else {
      m_mc[[m]][[1]] <- quote(lme4::glmer)
      m_mc[[m]]$control <- get_control_args(glmer = TRUE)               
    }
    y_mod[[m]] <- eval(m_mc[[m]], parent.frame())      
    
    # Error check: time_var is one of the longitudinal covariates
    # NB REMOVED, since difficult to check for time_var if it is 
    #   transformed (for example used in spline terms)
    #fm_nobars <- lme4::subbars(formula(y_mod[[m]]))
    #if (!time_var %in% rownames(attr(terms(fm_nobars), "factors")))
    #  stop(paste0("Variable '", time_var, "' does not appear in the ",
    #              "regression equation for the longitudinal submodel"))
    
    # Indicator of real or integer response vector
    y_is_real[m] <- check_response_real(family[[m]]$family)
    
    # Indicator of whether model includes a dispersion term
    y_has_dispersion[m] <- check_for_dispersion(family[[m]]$family)
    
    # Response vector and design matrix
    y[[m]] <- as.vector(lme4::getME(y_mod[[m]], "y"))
    x[[m]] <- as.matrix(lme4::getME(y_mod[[m]], "X"))
    y_has_intercept[m] <- grepl("(Intercept", colnames(x[[m]])[1L], fixed = TRUE)
    if (y_has_intercept[m]) {
      check_int <- check_intercept(family[[m]]$family, family[[m]]$link)
      y_has_intercept_unbound[m] <- check_int$unbound
      y_has_intercept_lobound[m] <- check_int$lobound
      y_has_intercept_upbound[m] <- check_int$upbound
    } else {
      y_has_intercept_unbound[m] <- 0L
      y_has_intercept_lobound[m] <- 0L
      y_has_intercept_upbound[m] <- 0L
    } 
    xtemp[[m]] <- if (y_has_intercept[m]) x[[m]][, -1L, drop=FALSE] else x[[m]]
    
    if (is.binomial(family[[m]]$family)) {
      if (NCOL(y[[m]]) == 1L) {
        if (is.numeric(y[[m]]) || is.logical(y[[m]])) 
          y[[m]] <- as.integer(y[[m]])
        if (is.factor(y[[m]])) 
          y[[m]] <- fac2bin(y[[m]])
        if (!all(y[[m]] %in% c(0L, 1L))) 
          stop("y values must be 0 or 1 for bernoulli regression.")
        trials[[m]] <- rep(0L, length(y[[m]]))
      } else {
        if (!isTRUE(NCOL(y[[m]]) == 2L))
          stop("y should either be a vector or a matrix 1 or 2 columns.")
        trials[[m]] <- as.integer(y[[m]][, 1L] + y[[m]][, 2L])
        y[[m]] <- as.integer(y[[m]][, 1L])
      }
    } else trials[[m]] <- rep(0L, length(y[[m]]))
    
    # Random effect terms
    Ztlist[[m]] <- lme4::getME(y_mod[[m]], "Ztlist")
    y_cnms[[m]]  <- lme4::getME(y_mod[[m]], "cnms")
    y_flist[[m]] <- lme4::getME(y_mod[[m]], "flist")
    y_offset[[m]] <- lme4::getME(y_mod[[m]], "offset")
    
    # Centred design matrix, if required
    if (centreLong) {
      y_centre[m] <- 1L
      xbar[[m]] <- colMeans(xtemp[[m]])
      xtemp[[m]] <- sweep(xtemp[[m]], 2, xbar[[m]], FUN = "-")
    }
    
    # Reorder y, X, Z if bernoulli (zeros first)
    if (is.binomial(family[[m]]$family) && all(y[[m]] %in% 0:1)) {      
      ord[[m]] <- order(y[[m]])
      y[[m]] <- y[[m]][ord[[m]]]
      trials[[m]] <- trials[[m]][ord[[m]]]
      xtemp[[m]] <- xtemp[[m]][ord[[m]], , drop = FALSE]  
      Ztlist[[m]] <- lapply(Ztlist[[m]], function(x) x[, ord[[m]], drop = FALSE])
      y_N01[[m]] <- sapply(0:1, function(x) sum(y[[m]] == x))
    } else y_N01[[m]] <- rep(0L, 2)  # dud entry if not bernoulli
    
    # Dimensions
    y_N[m] <- NROW(xtemp[[m]])
    y_real_N[m] <- if (y_is_real[m]) y_N[m] else 0L
    y_int_N[m] <- if (!y_is_real[m]) y_N[m] else 0L
    y_K[m] <- NCOL(xtemp[[m]])
    
    # Update formula if using splines or other data dependent predictors
    y_vars[[m]] <- get_formvars(y_mod[[m]])
    
    # Model based initial values or informative priors
    if (y_has_intercept_unbound[m]) {
      y_beta[[m]] <- lme4::fixef(y_mod[[m]])[-1L]
      se_y_beta[[m]] <- summary(y_mod[[m]])$coefficients[, "Std. Error"][-1L]
      markun <- sum(y_has_intercept_unbound[1:m])
      y_gamma_unbound[[markun]] <- fixef(y_mod[[m]])[1L]
      if (centreLong) 
        y_gamma_unbound[[markun]] <- y_gamma_unbound[[markun]] - xbar[[m]] %*% y_beta[[m]]
    } else if (y_has_intercept_lobound[m]) {
      y_beta[[m]] <- lme4::fixef(y_mod[[m]])[-1L]
      se_y_beta[[m]] <- summary(y_mod[[m]])$coefficients[, "Std. Error"][-1L]
      marklo <- sum(y_has_intercept_lobound[1:m])        
      y_gamma_lobound[[marklo]] <- fixef(y_mod[[m]])[1L]
      if (centreLong) 
        y_gamma_lobound[[marklo]] <- y_gamma_lobound[[marklo]] - xbar[[m]] %*% y_beta[[m]]
    } else if (y_has_intercept_upbound[m]) {
      y_beta[[m]] <- lme4::fixef(y_mod[[m]])[-1L]
      se_y_beta[[m]] <- summary(y_mod[[m]])$coefficients[, "Std. Error"][-1L]
      markup <- sum(y_has_intercept_upbound[1:m])        
      y_gamma_upbound[[markup]] <- fixef(y_mod[[m]])[1L]
      if (centreLong) 
        y_gamma_upbound[[markup]] <- y_gamma_upbound[[markup]] - xbar[[m]] %*% y_beta[[m]]
    } else {
      y_beta[[m]] <- lme4::fixef(y_mod[[m]])
      se_y_beta[[m]] <- summary(y_mod[[m]])$coefficients[, "Std. Error"]
    }
    if (y_has_intercept[m]) {
      tmp_gamma <- if (!centreLong) fixef(y_mod[[m]])[1L] else
        fixef(y_mod[[m]])[1L] - xbar[[m]] %*% y_beta[[m]]
      y_gamma <- c(y_gamma, tmp_gamma)
      se_y_gamma <- c(se_y_gamma, summary(y_mod[[m]])$coefficients[, "Std. Error"][1L])
    }
    vc <- lme4::VarCorr(y_mod[[m]])[[1]]
    b_Cov[[m]] <- vc
    sd_b[[m]] <- attr(vc, "stddev")
    b_Corr[[m]] <- attr(vc, "correlation")
    if (y_has_dispersion[m]) {
      disp_mark <- sum(y_has_dispersion[1:m])
      y_dispersion[disp_mark] <- sigma(y_mod[[m]])
    }
    
  }
  
  # Sum dimensions across all longitudinal submodels
  sum_y_N <- sum(y_N)
  sum_y_real_N <- sum(y_real_N)
  sum_y_int_N <- sum(y_int_N)
  sum_y_K <- sum(y_K)
  sum_y_has_intercept <- sum(y_has_intercept)
  sum_y_has_intercept_unbound <- sum(y_has_intercept_unbound)
  sum_y_has_intercept_lobound <- sum(y_has_intercept_lobound)
  sum_y_has_intercept_upbound <- sum(y_has_intercept_upbound)
  sum_y_has_dispersion <- sum(y_has_dispersion)
  
  # Indexing for binded response vector, design matrix, weights, etc
  y_beg <- sapply(1:M, function(m) sum(y_N[0:(m-1)]) + 1)
  y_end <- sapply(1:M, function(m) sum(y_N[0:m]))  
  y_real_beg <- sapply(1:M, function(m) if (y_is_real[m]) sum(y_real_N[0:(m-1)]) + 1L else 0L)
  y_real_end <- sapply(1:M, function(m) if (y_is_real[m]) sum(y_real_N[0:m]) else 0L)  
  y_int_beg <- sapply(1:M, function(m) if (!y_is_real[m]) sum(y_int_N[0:(m-1)]) + 1L else 0L)
  y_int_end <- sapply(1:M, function(m) if (!y_is_real[m]) sum(y_int_N[0:m]) else 0L)  
  
  # Additional error checks
  id_var <- check_id_var(id_var, y_cnms)
  id_list <- check_id_list(id_var, y_flist)
  
  # Construct weights
  if (has_weights) {
    if ((!is.data.frame(weights)) || (!ncol(weights) == 2))
      stop("'weights' argument should be a data frame with two columns: the first ",
           "containing patient IDs, the second containing their corresponding ",
           "weights.", call. = FALSE)
    if (!id_var %in% colnames(weights))
      stop("The data frame supplied in the 'weights' argument should have a ",
           "column named ", id_var, call. = FALSE)
    weight_var <- setdiff(colnames(weights), id_var)
    
    wts <- weights[[weight_var]]
    if (!is.numeric(wts)) 
      stop("The weights supplied must be numeric.", call. = FALSE)
    if (any(wts < 0)) 
      stop("Negative weights are not allowed.", call. = FALSE)
    
    # Check only one weight per ID
    n_weights_per_id <- tapply(weights[[weight_var]], weights[[id_var]], length)
    if (!all(n_weights_per_id == 1L))
      stop("The data frame supplied in the 'weights' argument should only have ",
           "one row (ie, one weight) per patient ID.", call. = FALSE)
    
    # Check for IDs with no weight supplied
    sel <- which(!id_list %in% factor(weights[[id_var]]))
    if (length(sel)) {
      if (length(sel) > 30L) sel <- sel[1:30]
      stop(paste0("The following patient IDs are used in fitting the model, but ",
                  "do not have weights supplied via the 'weights' argument: ",
                  paste(id_list[sel], collapse = ", ")), call. = FALSE)
    }
    for (m in 1:M) {
      # Obtain length and ordering of weights vector using flist
      flist_df <- data.frame(id = y_flist[[m]][[id_var]])
      weights_df <- merge(flist_df, weights, 
                          by.x = "id", by.y = id_var, sort = FALSE)
      y_weights[[m]] <- weights_df[[weight_var]]    
      # Reorder weights if bernoulli
      if (is.binomial(family[[m]]$family) && all(y[[m]] %in% 0:1))      
        y_weights[[m]] <- y_weights[[m]][ord[[m]]]
    }
  } else y_weights <- lapply(1:M, function(m) rep(0.0, length(y[[m]])))
  
  # Construct single cnms list for all longitudinal submodels
  y_cnms_nms <- lapply(y_cnms, names)
  cnms_nms <- unique(unlist(y_cnms_nms))
  cnms <- lapply(seq_along(cnms_nms), function(i) {
    nm <- cnms_nms[i]
    unlist(lapply(1:length(y_cnms), function(m) 
      if (nm %in% y_cnms_nms[[m]]) paste0("Long", m, "|", y_cnms[[m]][[nm]])))
  })
  names(cnms) <- cnms_nms
  
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
  
  # Require intercept for certain family and link combinations
  lapply(1:M, function(x) {
    if (!y_has_intercept[x]) {
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
  
  is_lmer <- lapply(family, is.lmer)
  
  #=========================
  # Data for event submodel
  #=========================
  
  # Items to store for event submodel
  e_beta <- c()
  
  # Baseline hazard
  base_haz <- match.arg(base_haz)
  base_haz_weibull <- (base_haz == "weibull")
  base_haz_piecewise <- (base_haz == "piecewise")
  base_haz_bs <- (base_haz == "bs")
  
  if (!(base_haz_bs || base_haz_piecewise)) { # not bs or piecewise
    if (!missing(df)) {
      warning("'df' will be ignored since 'base_haz' was not set ",
              "to B-splines or piecewise constant.", 
              immediate. = TRUE, call. = FALSE)
    }
    if (!missing(knots)) {
      warning("'knots' will be ignored since 'base_haz' was not set ",
              "to B-splines or piecewise constant.", 
              immediate. = TRUE, call. = FALSE)
    }
    df <- knots <- bs_df <- piecewise_df <- NULL    
  } else {  # bs  or piecewise
    if ((!missing(df)) && (!missing(knots))) {
      # both specified
      stop("Cannot specify both 'df' and 'knots'.", call. = FALSE)
    } else if (missing(df) && missing(knots)) {
      # neither specified -- use default df
      df <- bs_df <- piecewise_df <- 6L
      knots <- NULL
    } else if ((!missing(df)) && (missing(knots))) {
      # only df specified
      if (base_haz_bs) {
        if (df < 3)
          stop("'df' must be at least 3 for B-splines.")
        df <- bs_df <- df + 1
      } else if (base_haz_piecewise) {
        piecewise_df <- df
      }
      knots <- NULL
    } else if ((!missing(knots)) && (missing(df))) {
      # only knots specified
      if (any(knots < 0))
        stop("'knots' must be non-negative.")
      df <- NULL
      bs_df <- length(knots) + 4
      piecewise_df <- length(knots) + 1
    } else stop("Bug found: unable to reconcile 'df' and ",
                "'knots' arguments.", call. = FALSE)
  }
  
  # Set up model frame for event submodel
  if (!id_var %in% colnames(dataEvent))
    stop(paste0("Variable '", id_var, "' must be appear in dataEvent"),
         call. = FALSE)
  e_mc[[1]] <- quote(survival::coxph) 
  e_mc$x <- TRUE
  e_mod <- eval(e_mc, parent.frame())
  e_fr <- e_mf <- expand.model.frame(e_mod, id_var, na.expand = TRUE)
  e_mf <- cbind(unclass(e_mf[,1]), e_mf[, -1, drop = FALSE])
  
  # Check ID sorting
  e_id_list <- factor(unique(e_mf[, id_var]))
  if (!identical(id_list, e_id_list))
    stop("'dataEvent' needs to be sorted by the subject ",
         "ID/grouping variable", call. = FALSE)
  
  e_y <- e_mod$y
  
  # Entry and exit times
  entrytime <- rep(0, length(id_list)) # entry times current assumed to be zero for all individuals
  if (attr(e_y, "type") == "counting") {
    tvc         <- TRUE
    mf_event    <- do.call(rbind, lapply(split(e_mf, e_mf[, id_var]), function(d) d[which.max(d[,"stop"]), ]))
    eventtime   <- mf_event[["stop"]]
  } else if (attr(e_y, "type") == "right") {
    tvc         <- FALSE 
    mf_event    <- e_mf
    eventtime   <- mf_event[["time"]]
  } else stop("Only 'right' or 'counting' type Surv objects are allowed 
               on the LHS of the event submodel formula")
  
  # Event indicator and ID list
  d           <- mf_event[["status"]]  
  flist_event <- mf_event[[id_var]]
  names(eventtime) <- names(d) <- flist_event
  
  # Unstandardised quadrature points
  quadpoint <- lapply(quadpoints$points, unstandardise_quadpoints, entrytime, eventtime)
  t_q <- c(list(eventtime), quadpoint)
  names(quadpoint) <- paste0("quadpoint", seq(quadnodes))
  names(t_q) <- c("eventtime", names(quadpoint))
  
  # Obtain design matrix at event times and unstandardised quadrature points
  
  if (tvc) {  # time varying covariates in event model
    
    # Model frame at event times
    e_mf           <- data.table(e_mf, key = c(id_var, "start"))
    e_mf[["start"]] <- as.numeric(e_mf[["start"]])
    e_mf_eventtime <- e_mf[, .SD[.N], by = get(id_var)]
    e_mf_eventtime <- e_mf_eventtime[, get := NULL]   
    
    # Model frame corresponding to observation times which are 
    #   as close as possible to the unstandardised quadrature points                      
    e_mf_q  <- do.call(rbind, lapply(quadpoint, FUN = function(x)
      e_mf[data.table::SJ(flist_event, x), roll = TRUE, rollends = c(TRUE, TRUE)]))
    
    # Model frame evaluated at both event times and quadrature points
    e_mf_q <- rbind(e_mf_eventtime, e_mf_q)
    
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are time varying covariates in the event submodel and
    #   therefore the design matrix differs depending on the quadrature point 
    fm_RHS <- delete.response(terms(e_mod))
    e_x_quadtime   <- model.matrix(fm_RHS, data = e_mf_q)
    
  } else {  # no time varying covariates in event model
    
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are no time varying covariates in the event submodel and
    #   therefore the design matrix is identical at event time and at all
    #   quadrature points
    e_x_quadtime   <- do.call(rbind, lapply(1:(quadnodes + 1), FUN = function(x) e_mod$x))
    
  }
  
  # Evaluate spline basis (knots, df, etc) based on distributionof observed event times
  if (base_haz_bs) 
    bs_basis <- splines::bs(eventtime[(d > 0)], df = df, knots = knots, 
                            Boundary.knots = c(0, max(eventtime)), intercept = TRUE)
  
  # Evaluate cut points for piecewise constant
  if (base_haz_piecewise) {
    if (is.null(knots)) {
      knots <- quantile(eventtime[(d > 0)], probs = seq(0, 1, 1 / df))
      knots[[1]] <- 0
      knots[[length(knots)]] <- max(eventtime)
    } else {
      if (!is.numeric(knots))
        stop("'knots' vector must be numeric", call. = FALSE)
      if (any(knots < 0))
        stop("'knots' must be positive", call. = FALSE)
      if (any(knots > max(eventtime)))
        stop("'knots' cannot be greater than the largest event ",
             "time", call. = FALSE)
      knots <- c(0, knots, max(eventtime))
    }
  }
  
  # Incorporate intercept term (since Cox model does not have intercept)
  # -- depends on baseline hazard
  if (base_haz_weibull & (!"(Intercept)" %in% colnames(e_x_quadtime)))
    e_x_quadtime  <- cbind("(Intercept)" = rep(1, NROW(e_x_quadtime)), e_x_quadtime)
  
  # Centering of design matrix for event model
  e_x <- as.matrix(e_x_quadtime) 
  e_has_intercept <- if (length(colnames(e_x))) 
    grepl("(Intercept", colnames(e_x)[1L], fixed = TRUE) else FALSE
  e_xtemp <- if (e_has_intercept) e_x[, -1L, drop=FALSE] else e_x
  if (centreEvent) {
    e_xbar <- colMeans(e_xtemp)
    e_xtemp <- sweep(e_xtemp, 2, e_xbar, FUN = "-")
  }
  
    
  # Some additional bits -- NB weights here are quadrature weights, not prior weights
  e_K <- NCOL(e_xtemp)
  Npat <- length(eventtime)
  quadweight <- unlist(lapply(quadpoints$weights, unstandardise_quadweights, entrytime, eventtime))
  
  # Construct prior weights for event submodel
  if (has_weights) {
    flist_df <- data.frame(id = flist_event)
    weights_df <- merge(flist_df, weights, 
                        by.x = "id", by.y = id_var)
    e_weights <- weights_df[[weight_var]]
    e_weights_rep <- rep(e_weights, times = quadnodes)
  } else {
    e_weights <- rep(0.0, Npat)
    e_weights_rep <- rep(0.0, Npat * quadnodes)
  }
  
  # Error checks for the ID variable
  if (!identical(id_list, factor(sort(unique(flist_event)))))
    stop("The patient IDs (levels of the grouping factor) included ",
         "in the longitudinal and event submodels do not match")
  if (!identical(length(id_list), length(eventtime)))
    stop("The number of patients differs between the longitudinal and ",
         "event submodels. Perhaps you intended to use 'start/stop' notation ",
         "for the Surv() object.")
  
  # Model based initial values
  e_beta <- e_mod$coef
  se_e_beta <- sqrt(diag(e_mod$var))
  
  #================================
  # Data for association structure
  #================================
  
  # Time shift used for numerically calculating derivative of linear predictor 
  # or expected value of longitudinal outcome using one-sided difference
  eps <- 1E-5
  
  # Check association structure
  supported_assocs <- c("null", 
                        "etavalue", "muvalue", 
                        "etaslope", "muslope", 
                        "etalag",   "mulag",
                        "etaauc",   "muauc",
                        "shared_b", "shared_coef",
                        "etavalue_interact", "muvalue_interact",
                        "etaslope_interact", "muslope_interact")
  assoc <- lapply(1:M, validate_assoc, 
                  assoc, y_cnms, supported_assocs, id_var, x)
  assoc <- list_nms(assoc, M)
  
  # Indicator of each association type, for each longitudinal submodel
  has_assoc <- sapply(supported_assocs, function(x) 
    sapply(assoc, function(y) as.integer(y[[x]])), simplify = FALSE)
  
  which_b_zindex <- lapply(assoc, `[[`, "which_b_zindex")
  which_coef_zindex <- lapply(assoc, `[[`, "which_coef_zindex")
  which_coef_xindex <- lapply(assoc, `[[`, "which_coef_xindex")
  size_which_b <- sapply(which_b_zindex, length)
  size_which_coef <- sapply(which_coef_zindex, length)
  a_K <- get_num_assoc_pars(has_assoc, which_b_zindex, which_coef_zindex)  # doesn't include parameters for interaction terms (added on later)
  
  # Unstandardised quadrature nodes for AUC association structure
  auc_quadnodes <- 15
  auc_quadpoints <- get_quadpoints(auc_quadnodes)
  auc_quadweight <- unlist(
    lapply(t_q, function(x) 
      lapply(x, function(y) 
        lapply(auc_quadpoints$weights, unstandardise_quadweights, 0, y))))
  
  #==================================================
  # Longitudinal submodel: calculate design matrices 
  # at the event times and quadrature points
  #==================================================
  
  # Items to store for each longitudinal submodel
  y_mod_q         <- list()   # fitted long. submodels at event times and quadrature points
  y_fr            <- list()   # model frame with the addition of time_var
  xq              <- list()   # design matrix before removing intercept and centering
  xqtemp          <- list()   # design matrix after removing intercept and possibly centred
  Ztlistq         <- list()   # random effects matrix at event and quad times
  y_cnmsq         <- list() 
  y_flistq        <- list()
  
  y_mod_q_eps     <- list()   # fitted long. submodels under time shift of epsilon
  xq_eps          <- list()   
  xqtemp_eps      <- list()  
  Ztlistq_eps     <- list()
  y_cnmsq_eps     <- list()
  y_flistq_eps    <- list()
  
  y_mod_q_lag     <- list()   # fitted long. submodels under time lag
  xq_lag          <- list()  
  xqtemp_lag      <- list()  
  Ztlistq_lag     <- list()  
  y_cnmsq_lag     <- list()
  y_flistq_lag    <- list()
  
  y_mod_q_auc     <- list()   # fitted long. submodels at subquadrature points
  xq_auc          <- list()  
  xqtemp_auc      <- list()  
  Ztlistq_auc     <- list()  
  y_cnmsq_auc     <- list()
  y_flistq_auc    <- list()  
  
  xq_int          <- list()   # design matrix for covariates to interact with etavalue, etaslope, etc
  xqtemp_int      <- list()
  a_K_int         <- list()
  
  # Set up a second longitudinal model frame which includes the time variable
  for (m in 1:M) {
    m_mc_temp         <- m_mc[[m]]
    m_mc_temp[[1]]    <- quote(lme4::glFormula)    
    m_mc_temp$formula <- do.call(update.formula, list(m_mc[[m]]$formula, paste0("~ . +", time_var))) 
    m_mc_temp$control <- get_control_args(glmer = !is_lmer[[m]], norank = TRUE)
    
    # Obtain model frame with time_var definitely included
    y_mod_wtime <- eval(m_mc_temp, parent.frame())      
    y_fr[[m]] <- y_mod_wtime$fr
    
    # Convert model frame to data.table
    mf <- data.table::data.table(y_mod_wtime$fr, key = c(id_var, time_var))
    mf[[time_var]] <- as.numeric(mf[[time_var]])
    
    # Identify row in longitudinal data closest to event time or quadrature point
    #   NB if the quadrature point is earlier than the first observation time, 
    #   then covariates values are carried back to avoid missing values.
    #   In any other case, the 
    #   observed covariates values from the most recent observation time
    #   preceeding the quadrature point are carried forward to represent the 
    #   covariate value(s) at the quadrature point. (To avoid missingness  
    #   there is no limit on how far forwards or how far backwards covariate 
    #   values can be carried). If no time varying covariates are present in
    #   the longitudinal submodel (other than the time variable) then nothing 
    #   is carried forward or backward.
    mf_q <- do.call(rbind, lapply(t_q, FUN = function(x) 
      mf[data.table::SJ(flist_event, x), roll = TRUE, rollends = c(TRUE, TRUE)]))
    
    # Obtain long design matrix evaluated at event times and quadrature points
    names(mf_q)[names(mf_q) == "t_q"] <- time_var
    m_mc_temp$formula <- m_mc[[m]]$formula  # return to original formula
    m_mc_temp$formula <- use_these_vars(y_mod[[m]], 
                                        c(y_vars[[m]]$formvars$fixed[-1], y_vars[[m]]$formvars$random[-1]), 
                                        c(y_vars[[m]]$predvars$fixed[-1], y_vars[[m]]$predvars$random[-1]))
    m_mc_temp$control <- m_mc[[m]]$control  # return to original control args
    m_mc_temp$data    <- mf_q               # data at event and quadrature times
    y_mod_q[[m]] <- eval(m_mc_temp, parent.frame())     
    
    xq[[m]] <- as.matrix(y_mod_q[[m]]$X)
    xqtemp[[m]] <- if (y_has_intercept[m]) xq[[m]][, -1L, drop=FALSE] else xq[[m]]  
    Ztlistq[[m]] <- y_mod_q[[m]]$reTrms$Ztlist
    y_cnmsq[[m]] <- y_mod_q[[m]]$reTrms$cnms
    y_flistq[[m]] <- y_mod_q[[m]]$reTrms$flist
    
    #Needs working out to appropriately deal with offsets??
    #offset_quadtime <- model.offset(mod_quadtime$fr) %ORifNULL% double(0)
    
    # Centering of design matrix for longitudinal model at event times
    # and quadrature times 
    if (centreLong) xqtemp[[m]] <- sweep(xqtemp[[m]], 2, xbar[[m]], FUN = "-")
    
    # If association structure is based on slope, then calculate design 
    # matrices under a time shift of epsilon
    if (sum(has_assoc$etaslope, has_assoc$muslope, has_assoc$etaslope_interact, has_assoc$muslope_interact)) {
      mf_q_eps <- mf_q
      mf_q_eps[[time_var]] <- mf_q_eps[[time_var]] + eps
      m_mc_temp_eps <- m_mc_temp
      m_mc_temp_eps$data <- mf_q_eps
      y_mod_q_eps[[m]] <- eval(m_mc_temp_eps, parent.frame())
      xq_eps[[m]] <- as.matrix(y_mod_q_eps[[m]]$X)
      xqtemp_eps[[m]] <- if (y_has_intercept[m]) xq_eps[[m]][, -1L, drop=FALSE] else xq_eps[[m]]  
      if (centreLong) xqtemp_eps[[m]] <- sweep(xqtemp_eps[[m]], 2, xbar[[m]], FUN = "-")
      Ztlistq_eps[[m]] <- y_mod_q_eps[[m]]$reTrms$Ztlist
      y_cnmsq_eps[[m]] <- y_mod_q_eps[[m]]$reTrms$cnms
      y_flistq_eps[[m]] <- y_mod_q_eps[[m]]$reTrms$flist
    } 
    
    # If association structure is based on a time lag, then calculate design 
    # matrices under the specified time lag
    if (sum(has_assoc$etalag, has_assoc$mulag)) {
      t_q_lag <- lapply(t_q, function(x) {
        tmp <- x - assoc[[m]]$which_lag
        tmp[tmp < 0] <- 0.0  # use baseline where lagged t is before baseline
        tmp
      })
      mf_q_lag <- do.call(rbind, lapply(t_q_lag, FUN = function(x) 
        mf[data.table::SJ(flist_event, x), roll = TRUE, rollends = c(TRUE, TRUE)]))
      m_mc_temp_lag <- m_mc_temp
      m_mc_temp_lag$data <- mf_q_lag
      y_mod_q_lag[[m]] <- eval(m_mc_temp_lag, parent.frame())
      xq_lag[[m]] <- as.matrix(y_mod_q_lag[[m]]$X)
      xqtemp_lag[[m]] <- if (y_has_intercept[m]) xq_lag[[m]][, -1L, drop=FALSE] else xq_lag[[m]]  
      if (centreLong) xqtemp_lag[[m]] <- sweep(xqtemp_lag[[m]], 2, xbar[[m]], FUN = "-")
      Ztlistq_lag[[m]] <- y_mod_q_lag[[m]]$reTrms$Ztlist
      y_cnmsq_lag[[m]] <- y_mod_q_lag[[m]]$reTrms$cnms
      y_flistq_lag[[m]] <- y_mod_q_lag[[m]]$reTrms$flist
    } 
    
    # If association structure is based on area under the marker trajectory, then 
    # calculate design matrices at the subquadrature points
    if (sum(has_assoc$etaauc, has_assoc$muauc)) {
      # Return a design matrix that is (quadnodes * auc_quadnodes * Npat) rows
      t_q_auc <- lapply(t_q, function(x) 
        unlist(lapply(x, function(y) lapply(auc_quadpoints$points, unstandardise_quadpoints, 0, y))))
      mf_q_auc <- do.call(rbind, lapply(t_q_auc, function(x)
          mf[data.table::SJ(rep(flist_event, each = auc_quadnodes), x), roll = TRUE, rollends = c(TRUE, TRUE)]))
      m_mc_temp_auc <- m_mc_temp
      m_mc_temp_auc$data <- mf_q_auc
      y_mod_q_auc[[m]] <- eval(m_mc_temp_auc, parent.frame())
      xq_auc[[m]] <- as.matrix(y_mod_q_auc[[m]]$X)
      xqtemp_auc[[m]] <- if (y_has_intercept[m]) xq_auc[[m]][, -1L, drop=FALSE] else xq_auc[[m]]  
      if (centreLong) xqtemp_auc[[m]] <- sweep(xqtemp_auc[[m]], 2, xbar[[m]], FUN = "-")
      Ztlistq_auc[[m]] <- y_mod_q_auc[[m]]$reTrms$Ztlist
      y_cnmsq_auc[[m]] <- y_mod_q_auc[[m]]$reTrms$cnms
      y_flistq_auc[[m]] <- y_mod_q_auc[[m]]$reTrms$flist
    }      
    
    # If association structure is based on interactions with data, then calculate 
    # the design matrix which will be multiplied by etavalue, etaslope, muvalue or muslope
    xq_int[[m]] <- sapply(
      c("etavalue_interact", "etaslope_interact", 
        "muvalue_interact", "muslope_interact"),
      function(x) {
        if (has_assoc[[x]][[m]]) {
          fm <- assoc[[m]]$formula_interact[[x]]
          vars <- rownames(attr(terms.formula(fm), "factors"))
          if (is.null(vars))
            stop(paste0("No variables found in the formula specified for the '", x,
                        "' association structure.", call. = FALSE))
          ff <- ~ foo + bar
          gg <- parse(text = paste("~", paste(c(id_var, time_var), collapse = "+")))[[1L]]
          ff[[2L]][[2L]] <- fm[[2L]]
          ff[[2L]][[3L]] <- gg[[2L]]
          oldcall <- getCall(y_mod[[m]])
          naa <- if (is.null(interaction_data)) oldcall$na.action
          subset <- if (is.null(interaction_data)) oldcall$subset else NULL               
          df <- if (is.null(interaction_data)) eval(m_mc[[m]]$data) else interaction_data
          sel <- which(!vars %in% colnames(df))
          if (length(sel))
            stop(paste0("The following variables were specified in the formula for the '", x,
                        "' association structure, but they cannot be found in the data: ", 
                        paste0(vars, collapse = ", ")))
          mf <- eval(call("model.frame", ff, data = df, subset = subset, 
                          na.action = naa), envir = environment(formula(y_mod[[m]])))
          mf <- data.table::data.table(mf, key = c(id_var, time_var))
          mf[[time_var]] <- as.numeric(mf[[time_var]])
          mf_q_int <- do.call(rbind, lapply(t_q, FUN = function(x) 
            mf[data.table::SJ(flist_event, x), roll = TRUE, rollends = c(TRUE, TRUE)]))
          xq_int <- stats::model.matrix(fm, data = mf_q_int)
          if ("(Intercept)" %in% colnames(xq_int)) 
            xq_int <- xq_int[, -1L, drop = FALSE]
          if (!ncol(xq_int))
            stop(paste0("Bug found: A formula was specified for the '", x, "' association ", 
                        "structure, but the resulting design matrix has no columns."), call. = FALSE)
        } else xq_int <- matrix(0, length(unlist(t_q)), 0)
        xq_int
      }, simplify = FALSE, USE.NAMES = TRUE)
    a_K_int[[m]] <- sapply(xq_int[[m]], ncol)
    xqtemp_int[[m]] <- do.call(cbind, xq_int[[m]])
    
  }
  a_K <- a_K + sum(unlist(a_K_int))
  
  #=====================
  # Prior distributions
  #=====================
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus")
  ok_intercept_dists <- ok_dists[1:3]
  
  # Priors for longitudinal submodel(s)
  priorLong_scaled <- priorLong_ops$scaled
  priorLong_min_prior_scale <- priorLong_ops$min_prior_scale
  priorLong_scale_for_dispersion <- 
    as.array(maybe_broadcast(priorLong_ops$prior_scale_for_dispersion, sum_y_has_dispersion))
  
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
      # !!! Need to change this so that link can be specific to each submodel
      priorLong_scale <- 
        set_prior_scale(priorLong_scale, default = 2.5, 
                        link = family[[1]]$link)
    } else {
      priorLong_dist <- ifelse(priorLong_dist == "hs", 3L, 4L)
    }
    
    priorLong_df <- maybe_broadcast(priorLong_df, sum_y_K)
    priorLong_df <- as.array(pmin(.Machine$double.xmax, priorLong_df))
    priorLong_mean <- maybe_broadcast(priorLong_mean, sum_y_K)
    priorLong_mean <- as.array(priorLong_mean)
    priorLong_scale <- maybe_broadcast(priorLong_scale, sum_y_K)
    priorLong_scale <- as.array(priorLong_scale)
  }
  
  if (is.null(priorLong_intercept)) {
    priorLong_informative_for_intercept <- FALSE
    priorLong_dist_for_intercept <- 0L
    priorLong_mean_for_intercept <- as.array(rep(0, sum_y_has_intercept))
    priorLong_scale_for_intercept <- priorLong_df_for_intercept <- as.array(rep(1, sum_y_has_intercept))
  } else {
    if (!is.list(priorLong_intercept)) 
      stop("'priorLong_intercept' should be a named list.")
    priorLong_dist_for_intercept <- priorLong_intercept$dist
    priorLong_scale_for_intercept <- priorLong_intercept$scale
    priorLong_mean_for_intercept <- priorLong_intercept$location
    priorLong_df_for_intercept <- priorLong_intercept$df 
    priorLong_df_for_intercept[is.na(priorLong_df_for_intercept)] <- 1
    
    if (!priorLong_dist_for_intercept %in% unlist(ok_intercept_dists)) {
      stop("The prior distribution for the intercept should be one of ",
           paste(names(ok_intercept_dists), collapse = ", "))
    } else {
      priorLong_dist_for_intercept <- 
        ifelse(priorLong_dist_for_intercept == "normal", 1L, 2L)
      priorLong_scale_for_intercept <- 
        set_prior_scale(priorLong_scale_for_intercept, default = 10, 
                        link = family[[1]]$link)      
    }
    
    priorLong_df_for_intercept <- maybe_broadcast(priorLong_df_for_intercept, sum_y_has_intercept)
    priorLong_df_for_intercept <- as.array(pmin(.Machine$double.xmax, priorLong_df_for_intercept))
    priorLong_mean_for_intercept <- maybe_broadcast(priorLong_mean_for_intercept, sum_y_has_intercept)
    priorLong_mean_for_intercept <- as.array(priorLong_mean_for_intercept)
    priorLong_scale_for_intercept <- maybe_broadcast(priorLong_scale_for_intercept, sum_y_has_intercept)
    priorLong_scale_for_intercept <- as.array(priorLong_scale_for_intercept)
  }
  
  # Priors for event submodel
  priorEvent_scaled <- priorEvent_ops$scaled
  priorEvent_min_prior_scale <- priorEvent_ops$min_prior_scale
  priorEvent_scale_for_weibull <- priorEvent_ops$prior_scale_for_basehaz
  priorEvent_scale_for_bs <- priorEvent_ops$prior_scale_for_basehaz
  priorEvent_scale_for_piecewise <- priorEvent_ops$prior_scale_for_basehaz
  if (base_haz_bs) priorEvent_scale_for_bs <- 
    maybe_broadcast(priorEvent_scale_for_bs, bs_df)  
  if (base_haz_piecewise) priorEvent_scale_for_piecewise <- 
    maybe_broadcast(priorEvent_scale_for_piecewise, piecewise_df)  
  
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
      # !!! Need to think about whether 2 is appropriate default value here
      priorEvent_scale <- set_prior_scale(priorEvent_scale, default = 2, 
                                          link = "none")
    } else {
      priorEvent_dist <- ifelse(priorEvent_dist == "hs", 3L, 4L)
    }
    
    priorEvent_df <- maybe_broadcast(priorEvent_df, e_K)
    priorEvent_df <- as.array(pmin(.Machine$double.xmax, priorEvent_df))
    priorEvent_mean <- maybe_broadcast(priorEvent_mean, e_K)
    priorEvent_mean <- as.array(priorEvent_mean)
    priorEvent_scale <- maybe_broadcast(priorEvent_scale, e_K)
    priorEvent_scale <- as.array(priorEvent_scale)
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
    # !!! Need to think about whether 50 is appropriate default value here      
    priorEvent_scale_for_intercept <- 
      set_prior_scale(priorEvent_scale_for_intercept, default = 50, 
                      link = "none")
    priorEvent_df_for_intercept <- min(.Machine$double.xmax, priorEvent_df_for_intercept)
  }
  
  # Priors for association parameters
  priorAssoc_scaled <- priorAssoc_ops$scaled  # not currently used
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
      # !!! Need to potentially have appropriate default value here depending on 
      #     type of association structure
      priorAssoc_scale <- set_prior_scale(priorAssoc_scale, default = 25, 
                                          link = "none")      
    } else {
      priorAssoc_dist <- ifelse(priorAssoc_dist == "hs", 3L, 4L)
    }
    
    priorAssoc_df <- maybe_broadcast(priorAssoc_df, a_K)
    priorAssoc_df <- as.array(pmin(.Machine$double.xmax, priorAssoc_df))
    priorAssoc_mean <- maybe_broadcast(priorAssoc_mean, a_K)
    priorAssoc_mean <- as.array(priorAssoc_mean)
    priorAssoc_scale <- maybe_broadcast(priorAssoc_scale, a_K)
    priorAssoc_scale <- as.array(priorAssoc_scale)
  }
  
  # Minimum scaling of priors for longitudinal submodel(s)
  if (priorLong_scaled && priorLong_dist > 0L) {
    for (m in 1:M) {
      if (y_K[m] > 0L) {
        if (m == 1L) {
          mark_start <- 1 
          mark_end <- y_K[1]
        } else {
          mark_start <- sum(y_K[1:(m-1)]) + 1
          mark_end <- sum(y_K[1:m])
        }
        if (is_gaussian[[m]]) {
          ss <- 2 * sd(y[[m]])
          priorLong_scale[mark_start:mark_end] <- ss * priorLong_scale[mark_start:mark_end]
          priorLong_scale_for_intercept[[m]] <-  ss * priorLong_scale_for_intercept[[m]]
        }
        if (!QR) 
          priorLong_scale[mark_start:mark_end] <- 
            pmax(priorLong_min_prior_scale, priorLong_scale[mark_start:mark_end] / 
                   apply(xtemp[[m]], 2L, FUN = function(x) {
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
    as.array(pmin(.Machine$double.xmax, priorLong_scale_for_intercept))
  
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
    stop("QR decomposition is not yet supported by stan_jm")
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
    M = as.integer(M),
    Npat = as.integer(Npat),
    y_N = as.array(y_N), 
    y_real_N = as.array(y_real_N), 
    y_int_N = as.array(y_int_N), 
    y_N01 = as.array(t(sapply(y_N01, cbind))), 
    y_K = as.array(y_K), 
    sum_y_N = as.integer(sum_y_N),
    sum_y_real_N = as.integer(sum_y_real_N),
    sum_y_int_N = as.integer(sum_y_int_N),
    sum_y_K = as.integer(sum_y_K),
    e_K = as.integer(e_K),
    a_K = as.integer(a_K),
    quadnodes = as.integer(quadnodes),
    Npat_times_quadnodes = as.integer(Npat * quadnodes),
    sum_y_has_intercept = as.integer(sum_y_has_intercept), 
    sum_y_has_intercept_unbound = as.integer(sum_y_has_intercept_unbound), 
    sum_y_has_intercept_lobound = as.integer(sum_y_has_intercept_lobound), 
    sum_y_has_intercept_upbound = as.integer(sum_y_has_intercept_upbound), 
    sum_y_has_dispersion = as.integer(sum_y_has_dispersion),
    has_weights = as.integer(has_weights),
    
    # data for longitudinal submodel(s)
    link = as.array(link),
    y_centre = as.integer(centreLong),
    y_has_intercept = as.array(as.integer(y_has_intercept)),
    y_has_intercept_unbound = as.array(y_has_intercept_unbound),
    y_has_intercept_lobound = as.array(y_has_intercept_lobound),
    y_has_intercept_upbound = as.array(y_has_intercept_upbound),
    y_has_dispersion = as.array(as.numeric(y_has_dispersion)),
    y_real = as.array(as.numeric(unlist(y[y_is_real]))),
    y_int = as.array(as.integer(unlist(y[!y_is_real]))),
    y_beg = as.array(y_beg),  # indexing for combined response vector
    y_end = as.array(y_end),
    y_real_beg = as.array(y_real_beg),
    y_real_end = as.array(y_real_end),
    y_int_beg = as.array(y_int_beg), 
    y_int_end = as.array(y_int_end),
    y_weights = as.array(do.call(c, y_weights)),
    trials = as.array(as.integer(do.call(c, trials))),
    y_xbar = if (centreLong) as.array(do.call(c, xbar)) else double(0),
    y_X = as.array(as.matrix(Matrix::bdiag(xtemp))),
    
    # data for event submodel
    basehaz_weibull = as.integer(base_haz_weibull),
    basehaz_piecewise = as.integer(base_haz_piecewise),
    basehaz_bs = as.integer(base_haz_bs),
    bs_df = if (base_haz_bs) as.integer(bs_df) else 0,
    piecewise_df = if (base_haz_piecewise) as.integer(piecewise_df) else 0,
    e_centre = as.integer(centreEvent),
    e_has_intercept = as.integer(e_has_intercept),
    nrow_y_Xq = NROW(xqtemp[[1]]),
    nrow_e_Xq = NROW(e_xtemp),
    y_Xq = as.array(as.matrix(Matrix::bdiag(xqtemp))),
    e_Xq = e_xtemp,
    e_times = c(eventtime, unlist(quadpoint)),
    e_d = c(d, rep(1, length(unlist(quadpoint)))),
    e_xbar = if (centreEvent) as.array(e_xbar) else double(0),
    e_weights = as.array(e_weights),
    e_weights_rep = as.array(e_weights_rep),
    quadweight = as.array(quadweight),
    
    # data for association structure
    assoc = as.integer(a_K > 0L),
    has_assoc_ev = as.array(as.integer(has_assoc$etavalue)),
    has_assoc_es = as.array(as.integer(has_assoc$etaslope)),
    has_assoc_el = as.array(as.integer(has_assoc$etalag)),
    has_assoc_ea = as.array(as.integer(has_assoc$etaauc)),
    has_assoc_mv = as.array(as.integer(has_assoc$muvalue)),
    has_assoc_ms = as.array(as.integer(has_assoc$muslope)),
    has_assoc_ml = as.array(as.integer(has_assoc$mulag)),
    has_assoc_ma = as.array(as.integer(has_assoc$muauc)),
    has_assoc_evi = as.array(as.integer(has_assoc$etavalue_interact)),
    has_assoc_esi = as.array(as.integer(has_assoc$etaslope_interact)),
    has_assoc_mvi = as.array(as.integer(has_assoc$muvalue_interact)),
    has_assoc_msi = as.array(as.integer(has_assoc$muslope_interact)),    
    sum_has_assoc_ev = as.integer(sum(has_assoc$etavalue)),
    sum_has_assoc_es = as.integer(sum(has_assoc$etaslope)),
    sum_has_assoc_el = as.integer(sum(has_assoc$etalag)),
    sum_has_assoc_ea = as.integer(sum(has_assoc$etaauc)),
    sum_has_assoc_mv = as.integer(sum(has_assoc$muvalue)),
    sum_has_assoc_ms = as.integer(sum(has_assoc$muslope)),
    sum_has_assoc_ml = as.integer(sum(has_assoc$mulag)),
    sum_has_assoc_ma = as.integer(sum(has_assoc$muauc)),
    sum_has_assoc_evi = as.integer(sum(has_assoc$etavalue_interact)),
    sum_has_assoc_esi = as.integer(sum(has_assoc$etaslope_interact)),
    sum_has_assoc_mvi = as.integer(sum(has_assoc$muvalue_interact)),
    sum_has_assoc_msi = as.integer(sum(has_assoc$muslope_interact)),
    auc_quadnodes = as.integer(auc_quadnodes),
    auc_quadweight = as.array(auc_quadweight),
    nrow_y_Xq_auc = as.integer(auc_quadnodes * NROW(xqtemp[[1]])),
    Npat_times_auc_quadnodes = as.integer(Npat * auc_quadnodes),
    sum_size_which_b = as.integer(sum(size_which_b)),
    size_which_b = as.array(size_which_b),
    which_b_zindex = as.array(unlist(which_b_zindex)),
    sum_size_which_coef = as.integer(sum(size_which_coef)),
    size_which_coef = as.array(size_which_coef),
    which_coef_zindex = as.array(unlist(which_coef_zindex)),
    which_coef_xindex = as.array(unlist(which_coef_xindex)),
    
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
    priorLong_scale_for_dispersion = as.array(priorLong_scale_for_dispersion),
    priorEvent_scale_for_weibull = 
      if (base_haz_weibull) as.array(priorEvent_scale_for_weibull) else as.array(double(0)),
    priorEvent_scale_for_bs = 
      if (base_haz_bs) as.array(priorEvent_scale_for_bs) else as.array(double(0)),
    priorEvent_scale_for_piecewise = 
      if (base_haz_piecewise) as.array(priorEvent_scale_for_piecewise) else as.array(double(0)),
    
    prior_PD = as.integer(prior_PD)
  )  
  
  # data for random effects
  group <- lapply(1:M, function(x) {
    pad_reTrms(Ztlist = Ztlist[[x]], 
               cnms = y_cnms[[x]], 
               flist = y_flist[[x]])})
  Z     <- lapply(1:M, function(x) group[[x]]$Z)
  y_cnms <- lapply(1:M, function(x) group[[x]]$cnms)
  y_flist_padded <- lapply(1:M, function(x) group[[x]]$flist)
  t <- length(cnms_nms) # num. of unique grouping factors
  t_i <- which(cnms_nms == id_var) # index of patient-level grouping factor
  pmat <- matrix(0, t, M)
  for (i in 1:t) {
    for (j in 1:M) {
      pmat[i,j] <- length(y_cnms[[j]][[cnms_nms[i]]])
    }
  }
  p <- rowSums(pmat)
  lmat <- matrix(0, t, M)
  l <- c()
  for (i in 1:t) {
    for (j in 1:M) {
      lmat[i,j] <- nlevels(y_flist_padded[[j]][[cnms_nms[i]]])
    }
    l[i] <- max(lmat[i,])
    if (!all(lmat[i,] %in% c(0, l[i])))
      stop("The number of factor levels for each of the grouping factors ",
           "must be the same in each of the longitudinal submodels")
  }
  qmat <- l * pmat
  q1 <- rowSums(qmat)
  q2 <- colSums(qmat)
  
  # Names of clustering variables
  group_nms <- lapply(y_cnms, names)
  # Names of random effects and random coefficients
  b_nms <- character()
  g_nms <- character() 
  for (m in 1:M) {
    for (i in seq_along(group_nms[[m]])) {
      # !!! if you change this change .pp_data_mer_z() as well
      nm <- group_nms[[m]][i]
      nms_i <- paste(y_cnms[[m]][[nm]], nm)
      if (length(nms_i) == 1) {
        b_nms <- c(b_nms, paste0("Long", m, "|", nms_i, ":", levels(y_flist_padded[[m]][[nm]])))
      } else {
        b_nms <- c(b_nms, c(t(sapply(nms_i, function(x) 
          paste0("Long", m, "|", x, ":", levels(y_flist_padded[[m]][[nm]]))))))
      }
      g_nms <- c(g_nms, paste0("Long", m, "|", nms_i)) 
    }
  }
  standata$t <- t
  standata$t_i <- as.integer(t_i)
  standata$pmat <- as.array(pmat)
  standata$p <- as.array(p)
  standata$l <- as.array(l)
  standata$qmat <- as.array(qmat)
  standata$q1 <- as.array(q1)
  standata$q2 <- as.array(q2)
  standata$len_theta_L <- sum(choose(p, 2), p)
  Zmerge <- Matrix::bdiag(Z)
  standata$len_b <- as.integer(ncol(Zmerge))
  parts <- rstan::extract_sparse_parts(Zmerge)
  standata$num_non_zero <- as.integer(length(parts$w))
  standata$w <- parts$w
  standata$v <- parts$v
  standata$u <- as.array(parts$u)
  
  # data for random effects in GK quadrature
  groupq <- lapply(1:M, function(x) {
    pad_reTrms(Ztlist = Ztlistq[[x]], 
               cnms = y_cnmsq[[x]], 
               flist = y_flistq[[x]])})
  Zq <- lapply(1:M, function(x) groupq[[x]]$Z)
  Zqmerge <- Matrix::bdiag(Zq)
  parts_Zq <- rstan::extract_sparse_parts(Zqmerge)
  standata$num_non_zero_Zq <- as.integer(length(parts_Zq$w))
  standata$w_Zq <- parts_Zq$w
  standata$v_Zq <- parts_Zq$v
  standata$u_Zq <- as.array(parts_Zq$u)
  
  # data for calculating eta slope in GK quadrature 
  standata$eps <- eps  # time shift for numerically calculating derivative
  standata$y_Xq_eps <- if (sum(has_assoc$etaslope, has_assoc$muslope))
    as.array(as.matrix(Matrix::bdiag(xqtemp_eps))) else as.array(matrix(0,0,sum_y_K))
  if (length(Ztlistq_eps)) {
    groupq_eps <- lapply(1:M, function(x) {
      pad_reTrms(Ztlist = Ztlistq_eps[[x]], 
                 cnms = y_cnmsq_eps[[x]], 
                 flist = y_flistq_eps[[x]])})
    Zq_eps <- lapply(1:M, function(x) groupq_eps[[x]]$Z)
    Zq_eps_merge <- Matrix::bdiag(Zq_eps) 
  } else Zq_eps_merge <- matrix(0,0,0)
  parts_Zq_eps <- rstan::extract_sparse_parts(Zq_eps_merge)
  standata$num_non_zero_Zq_eps <- as.integer(length(parts_Zq_eps$w))
  standata$w_Zq_eps <- parts_Zq_eps$w
  standata$v_Zq_eps <- parts_Zq_eps$v
  standata$u_Zq_eps <- as.array(parts_Zq_eps$u)    
  
  # data for calculating eta lag in GK quadrature 
  standata$y_Xq_lag <- if (sum(has_assoc$etalag, has_assoc$mulag))
    as.array(as.matrix(Matrix::bdiag(xqtemp_lag))) else as.array(matrix(0,0,sum_y_K))
  if (length(Ztlistq_lag)) {
    groupq_lag <- lapply(1:M, function(x) {
      pad_reTrms(Ztlist = Ztlistq_lag[[x]], 
                 cnms = y_cnmsq_lag[[x]], 
                 flist = y_flistq_lag[[x]])})
    Zq_lag <- lapply(1:M, function(x) groupq_lag[[x]]$Z)
    Zq_lag_merge <- Matrix::bdiag(Zq_lag) 
  } else Zq_lag_merge <- matrix(0,0,0)
  parts_Zq_lag <- rstan::extract_sparse_parts(Zq_lag_merge)
  standata$num_non_zero_Zq_lag <- as.integer(length(parts_Zq_lag$w))
  standata$w_Zq_lag <- parts_Zq_lag$w
  standata$v_Zq_lag <- parts_Zq_lag$v
  standata$u_Zq_lag <- as.array(parts_Zq_lag$u)    
  
  # data for calculating eta auc in GK quadrature 
  standata$y_Xq_auc <- if (sum(has_assoc$etaauc, has_assoc$muauc))
    as.array(as.matrix(Matrix::bdiag(xqtemp_auc))) else as.array(matrix(0,0,sum_y_K))
  if (length(Ztlistq_auc)) {
    groupq_auc <- lapply(1:M, function(x) {
      pad_reTrms(Ztlist = Ztlistq_auc[[x]], 
                 cnms = y_cnmsq_auc[[x]], 
                 flist = y_flistq_auc[[x]])})
    Zq_auc <- lapply(1:M, function(x) groupq_auc[[x]]$Z)
    Zq_auc_merge <- Matrix::bdiag(Zq_auc) 
  } else Zq_auc_merge <- matrix(0,0,0)
  parts_Zq_auc <- rstan::extract_sparse_parts(Zq_auc_merge)
  standata$num_non_zero_Zq_auc <- as.integer(length(parts_Zq_auc$w))
  standata$w_Zq_auc <- parts_Zq_auc$w
  standata$v_Zq_auc <- parts_Zq_auc$v
  standata$u_Zq_auc <- as.array(parts_Zq_auc$u)    
  
  # data for calculating interactions in GK quadrature 
  standata$y_Xq_int <- as.array(t(as.matrix(do.call("cbind", (xqtemp_int)))))
  standata$a_K_int <- as.array(unlist(a_K_int))  # number of columns in xqtemp_int corresponding to each interaction type (etavalue, etaslope, muvalue, muslope) for each submodel
  standata$sum_a_K_int <-  as.integer(sum(unlist(a_K_int)))
  
  # hyperparameters for random effects model
  if (prior_covariance$dist == "decov") {
    decov_args <- prior_covariance
    standata$shape <- as.array(maybe_broadcast(decov_args$shape, t))
    standata$scale <- as.array(maybe_broadcast(decov_args$scale, t))
    standata$len_concentration <- sum(p[p > 1])
    standata$concentration <- 
      as.array(maybe_broadcast(decov_args$concentration, sum(p[p > 1])))
    standata$len_regularization <- sum(p > 1)
    standata$regularization <- 
      as.array(maybe_broadcast(decov_args$regularization, sum(p > 1))) 
  }
  
  standata$family <- as.array(sapply(1:M, function(x) {
    return_fam <- switch(family[[x]]$family, 
                         gaussian = 1L, 
                         Gamma = 2L,
                         inverse.gaussian = 3L,
                         binomial = 5L,
                         poisson = 6L,
                         "neg_binomial_2" = 7L)
    if (is_bernoulli[[x]]) return_fam <- 4L
    return_fam}))
  standata$any_fam_3 <- as.integer(any(standata$family == 3L))
  
  # B-splines baseline hazard
  standata$e_ns_times <- if (base_haz_bs) 
    as.array(predict(bs_basis, standata$e_times)) else 
      as.array(matrix(0,0,0))
  
  # Piecewise constant baseline hazard
  if (base_haz_piecewise) {
    e_times_quantiles <- cut(standata$e_times, knots, 
                             include.lowest = TRUE, labels = FALSE)
    tmp <- matrix(NA, length(e_times_quantiles), piecewise_df)
    for (i in 1:piecewise_df) tmp[, i] <- ifelse(e_times_quantiles == i, 1, 0)
  }
  standata$e_times_piecedummy <- if (base_haz_piecewise)
    as.array(tmp) else as.array(matrix(0,0,0))    
  
  
  #================
  # Initial values
  #================
  
  if (init == "model_based") {
    y_gamma_unbound <- unlist(y_gamma_unbound)
    y_gamma_lobound <- unlist(y_gamma_lobound)
    y_gamma_upbound <- unlist(y_gamma_upbound)
    y_z_beta <- (unlist(y_beta) - priorLong_mean) / priorLong_scale
    y_dispersion_unscaled <- y_dispersion / priorLong_scale_for_dispersion
    e_z_beta <- (e_beta - priorEvent_mean) / priorEvent_scale 
    
    y_hs <- if (priorLong_dist <= 2L) 0 
    else if (priorLong_dist == 3L) 2
    else if (priorLong_dist == 4L) 4
    e_hs <- if (priorEvent_dist <= 2L) 0 
    else if (priorEvent_dist == 3L) 2
    else if (priorEvent_dist == 4L) 4
    a_hs <- if (priorAssoc_dist <= 2L) 0 
    else if (priorAssoc_dist == 3L) 2
    else if (priorAssoc_dist == 4L) 4
    
    if (prior_covariance$dist == "decov") {
      len_z_T <- 0
      for (i in 1:t) {
        if (p[i] > 2) 
          for (j in 3:p[i]) len_z_T <- len_z_T + p[i] - 1;
      }
      # Not yet calculating initial values for parameters which
      # are part of the decomposed theta_L matrix
      tau <- c()
      rho <- c()
      normalised_zetas <- c()
      rho_mark <- 1
      for (i in 1:t) {
        L_b_Cov <- t(chol(as.matrix(Matrix::bdiag(b_Cov))))
        diag_L_b_Cov <- diag(L_b_Cov)
        trace <- sum(diag_L_b_Cov)  # equal to variance of RE if only one RE
        normalised_zetas_tmp <- diag_L_b_Cov / trace  # equal to 1 if only one random effect
        #tau[i] <- std_dev / standata$scale[i]
        tau[i] <- (sqrt(trace / p[i])) / standata$scale[i]
        std_dev1 <- sqrt(L_b_Cov[1,1])
        if (p[i] > 1) {
          std_dev2 <- sqrt(L_b_Cov[2,2])
          T21 <- L_b_Cov[2,1] / std_dev2
          rho[rho_mark] <- (T21 + 1) / 2
          rho_mark <- rho_mark + 1
          normalised_zetas <- c(normalised_zetas, normalised_zetas_tmp)
        }
      }
    }
    
    model_based_inits <- Filter(function(x) (!is.null(x)), c(
      list(
        y_gamma_unbound = if (sum_y_has_intercept_unbound) as.array(y_gamma_unbound) else double(0),
        y_gamma_lobound = if (sum_y_has_intercept_lobound) as.array(y_gamma_lobound) else double(0),
        y_gamma_upbound = if (sum_y_has_intercept_upbound) as.array(y_gamma_upbound) else double(0),
        y_z_beta = if (sum_y_K) as.array(y_z_beta) else double(0),
        y_dispersion_unscaled = if (sum_y_has_dispersion) as.array(y_dispersion_unscaled) else double(0),
        e_gamma = if (e_has_intercept) as.array(0) else double(0),
        e_z_beta = if (e_K) as.array(e_z_beta) else double(0),
        weibull_shape_unscaled = if (base_haz_weibull) 
          as.array(runif(1, 0.5, 3) / priorEvent_scale_for_weibull) else double(0),
        bs_coefs_unscaled = if (base_haz_bs) as.array(rep(0, bs_df)) else double(0),
        piecewise_coefs_unscaled = if (base_haz_piecewise) as.array(rep(0, piecewise_df)) else double(0),
        a_z_beta = if (a_K) as.array(rep(0, a_K)) else double(0),
        z_b = as.array(runif(standata$len_b, -0.5, 0.5)),
        y_global = as.array(runif(y_hs)),
        y_local = if (y_hs) matrix(runif(y_hs * sum_y_K), nrow = y_hs, ncol = sum_y_K)
        else matrix(0,0,0),
        e_global = as.array(runif(e_hs)),
        e_local = if (e_hs) matrix(runif(e_hs * e_K), nrow = e_hs, ncol = e_K)
        else matrix(0,0,0),
        a_global = as.array(runif(a_hs)),
        a_local = if (a_hs) matrix(runif(a_hs * a_K), nrow = a_hs, ncol = a_K)
        else matrix(0,0,0)      
      ),
      if (prior_covariance$dist == "decov") list(
        z_T = as.array(rep(sqrt(1/len_z_T), len_z_T)),
        rho = if ((sum(p) - t) > 0) as.array(rep(1 / (sum(p) - t + 1), (sum(p) - t))) else double(0),
        zeta = if (!is.null(normalised_zetas)) as.array(normalised_zetas) else double(0),
        tau = as.array(tau)
      )
    ))
    init <- function() model_based_inits
  }
  
  
  #===========
  # Fit model
  #===========
    
  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$jm
  pars <- c(if (sum_y_has_intercept_unbound) "y_gamma_unbound",
            if (sum_y_has_intercept_lobound) "y_gamma_lobound",
            if (sum_y_has_intercept_upbound) "y_gamma_upbound",
            #if (sum_y_has_intercept_unbound) "y_alpha_unbound", 
            #if (sum_y_has_intercept_lobound) "y_alpha_lobound", 
            #if (sum_y_has_intercept_upbound) "y_alpha_upbound", 
            if (sum_y_K) "y_beta",
            if (e_has_intercept) "e_gamma",
            if (e_K) "e_beta",
            if (a_K) "a_beta",
            if (Npat) "b_by_model",
            "y_dispersion", 
            if (base_haz_weibull) "weibull_shape",
            if (base_haz_piecewise) "piecewise_coefs",
            if (base_haz_bs) "bs_coefs")
  
  cat(paste0(if (M == 1L) "Uni" else "Multi", 
             "variate joint model specified\n"))
  if (algorithm == "sampling") {
    cat("\nPlease note the warmup phase may be much slower than",
        "later iterations!\n")             
    sampling_args <- set_sampling_args_for_jm(
      object = stanfit, 
      user_dots = list(...), 
      user_adapt_delta = adapt_delta,
      user_max_treedepth = max_treedepth,
      sum_p = sum(standata$p),
      data = standata, 
      pars = pars, 
      init = init,
      show_messages = FALSE)
    stanfit <- do.call(sampling, sampling_args)
  } else {
    # meanfield or fullrank vb
    stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                         algorithm = algorithm, ...)    
  }
  
  #  if (QR) {  # not yet implemented for stan_jm
  #    thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, 
  #                      permuted = FALSE)
  #    betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
  #    end <- utils::tail(dim(betas), 1L)
  #    for (chain in 1:end) for (param in 1:nrow(betas)) {
  #      stanfit@sim$samples[[chain]][[has_intercept + param]] <-
  #        if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
  #    }
  #  }
  
  # Names for coefs from submodel(s)
  int_nms <- unlist(lapply(1:M, function(x) 
    if (y_has_intercept[x]) paste0("Long", x, "|(Intercept)")))
  y_nms   <- unlist(lapply(1:M, function(x) 
    if (ncol(xtemp[[x]])) paste0("Long", x, "|", colnames(xtemp[[x]]))))
  e_nms   <- if (ncol(e_x)) paste0("Event|", colnames(e_x))    
  
  # Names for vector of association parameters
  a_nms <- character()  
  for (m in 1:M) {
    if (has_assoc$etavalue[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etavalue"))
    if (has_assoc$etavalue_interact[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etavalue:", colnames(xq_int[[m]][["etavalue_interact"]])))
    if (has_assoc$etaslope[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etaslope"))
    if (has_assoc$etaslope_interact[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etaslope:", colnames(xq_int[[m]][["etaslope_interact"]])))    
    if (has_assoc$etalag[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etalag"))
    if (has_assoc$etaauc[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etaauc"))
    if (has_assoc$muvalue[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue"))
    if (has_assoc$muvalue_interact[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue:", colnames(xq_int[[m]][["muvalue_interact"]])))    
    if (has_assoc$muslope[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muslope"))
    if (has_assoc$muslope_interact[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muslope:", colnames(xq_int[[m]][["muslope_interact"]])))    
    if (has_assoc$mulag[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|mulag"))
    if (has_assoc$muauc[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muauc"))
  }
  if (sum(size_which_b)) {
    temp_g_nms <- lapply(1:M, FUN = function(m) {
      all_nms <- paste0(paste0("Long", m, "|b["), y_cnms[[m]][[id_var]], "]")
      all_nms[which_b_zindex[[m]]]})
    a_nms <- c(a_nms, paste0("Assoc|", unlist(temp_g_nms)))
  }
  if (sum(size_which_coef)) {
    temp_g_nms <- lapply(1:M, FUN = function(m) {
      all_nms <- paste0(paste0("Long", m, "|coef["), y_cnms[[m]][[id_var]], "]")
      all_nms[which_coef_zindex[[m]]]})
    a_nms <- c(a_nms, paste0("Assoc|", unlist(temp_g_nms)))
  }
  
  # Names for vector of dispersion parameters
  d_nms <- character()  
  for (m in 1:M) {
    if (is.gaussian(famname[[m]]))   d_nms <- c(d_nms, paste0("Long", m,"|sigma"))
    else if (is.gamma(famname[[m]])) d_nms <- c(d_nms, paste0("Long", m,"|shape"))
    else if (is.ig(famname[[m]]))    d_nms <- c(d_nms, paste0("Long", m,"|lambda"))
    else if (is.nb(famname[[m]]))    d_nms <- c(d_nms, paste0("Long", m,"|overdispersion"))
  }
  
  new_names <- c(int_nms,
                 y_nms,
                 e_nms,
                 a_nms,                   
                 if (length(group)) c(paste0("b[", b_nms, "]")),
                 d_nms,
                 if (base_haz_weibull) "Event|weibull-shape",               
                 if (base_haz_piecewise) paste0("Event|basehaz-coef", seq(piecewise_df)),               
                 if (base_haz_bs) paste0("Event|basehaz-coef", seq(bs_df)),
                 "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  
  n_grps <- l - 1
  names(n_grps) <- cnms_nms  # n_grps is num. of levels within each grouping factor
  names(p) <- cnms_nms       # p is num. of variables within each grouping factor
  
  # Attributes of baseline haz for passing to fitted object
  if (base_haz_bs) {
    attr <- attributes(bs_basis)[c("degree", "knots", "Boundary.knots", "intercept")] 
  } else if (base_haz_piecewise) {
    attr <- list(df = piecewise_df, knots = knots)
  } else attr <- NULL
  
  # Undo ordering of matrices if bernoulli
  for (m in 1:M) {
    if (is_bernoulli[[m]]) {
      y[[m]] <- y[[m]][order(ord[[m]])]
      trials[[m]] <- trials[[m]][order(ord[[m]])]
      y_weights[[m]] <- y_weights[[m]][order(ord[[m]])]
      xtemp[[m]] <- xtemp[[m]][order(ord[[m]]), , drop = FALSE]  
      Z[[m]] <- Z[[m]][order(ord[[m]]), , drop = FALSE]
    }
  }
  
  #colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula = c(formulaLong, formulaEvent), 
               id_var, time_var, offset = NULL, quadnodes,
               base_haz = list(type = base_haz, attr = attr),
               M, cnms, y_N, y_cnms, y_flist, Npat, n_grps, assoc = lapply(has_assoc, as.logical), 
               fr = c(y_fr, list(e_fr)),
               x = lapply(1:M, function(i) 
                 if (getRversion() < "3.2.0") Matrix::cBind(x[[i]], Z[[i]]) else cbind2(x[[i]], Z[[i]])),
               xq = lapply(1:M, function(i) 
                 if (getRversion() < "3.2.0") Matrix::cBind(xq[[i]], Zq[[i]]) else cbind2(xq[[i]], Zq[[i]])),                 
               xq_eps = lapply(1:M, function(i) 
                 if (has_assoc$etaslope[m] || has_assoc$muslope[m]) {
                   if (getRversion() < "3.2.0") 
                     Matrix::cBind(xq_eps[[i]], Zq_eps[[i]]) else cbind2(xq_eps[[i]], Zq_eps[[i]])
                 } else NULL),                          
               y = y, e_x, eventtime, d, quadpoints = quadpoint, 
               epsilon = if (sum(has_assoc$etaslope, has_assoc$muslope)) eps else NULL,
               standata, dataLong, dataEvent, call, terms = NULL, model = NULL,                          
               prior.info = get_prior_info(call, formals()),
               na.action, algorithm = "sampling", init, glmod = y_mod, coxmod = e_mod)
  out <- stanjm(fit)
  
  return(out)
}


# Set arguments for sampling for stan_jm
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# @param object The stanfit object to use for sampling.
# @param user_dots The contents of \code{...} from the user's call to
#   the \code{stan_jm} modeling function.
# @param user_adapt_delta The value for \code{adapt_delta} specified by the
#   user.
# @param user_max_treedepth The value for \code{max_treedepth} specified by the
#   user.
# @param sum_p The total number of random effects in the joint model. Should
#   likely be passed as sum(standata$p)
# @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#   \code{user_dots} (e.g. \code{pars}, \code{init}, etc.)
# @return A list of arguments to use for the \code{args} argument for 
#   \code{do.call(sampling, args)}.
set_sampling_args_for_jm <- function(object, user_dots = list(), 
                                     user_adapt_delta = NULL, 
                                     user_max_treedepth = NULL, 
                                     sum_p = NULL, ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  
  if (is.null(sum_p) || (sum_p <= 0))
    stop("Bug found: sum_p should specify the total number ",
         "of random effects in the joint model (used for ",
         "determining the default adapt_delta")
  
  default_adapt_delta <- if (sum_p > 2) 0.85 else 0.80
  default_max_treedepth <- 9L
  
  if (!is.null(user_adapt_delta))
    args$control$adapt_delta <- user_adapt_delta else 
      if (is.null(args$control$adapt_delta))
        args$control$adapt_delta <- default_adapt_delta
  
  if (!is.null(user_max_treedepth))
    args$control$max_treedepth <- user_max_treedepth else
      if (is.null(args$control$max_treedepth))
        args$control$max_treedepth <- default_max_treedepth
  
  if (!"iter" %in% unms) args$iter <- 1000
  if (!"chains" %in% unms) args$chains <- 1
  #if (!"cores" %in% unms) args$cores <- parallel::detectCores()
  if (!"refresh" %in% unms) args$refresh <- args$iter / 25
  if (!"save_warmup" %in% unms) args$save_warmup <- TRUE  
  
  return(args)
}  


# Check the id_var argument is valid and is included appropriately in the
# formulas for each of the longitudinal submodels
#
# @param id_var The character string that the user specified for the id_var
#   argument -- will have been set to NULL if the argument was missing.
# @param y_cnms A list of length M with the cnms for each longitudinal submodel
# @return Returns the character string corresponding to the appropriate id_var.
#   This will either be the user specified id_var argument or the only grouping
#   factor.
check_id_var <- function(id_var, y_cnms) {
  len_cnms <- sapply(y_cnms, length)
  if (any(len_cnms > 1L)) {  # more than one grouping factor
    if (is.null(id_var)) {
      stop("'id_var' must be specified when using more than one grouping factor",
           call. = FALSE)
    } else {
      lapply(y_cnms, function(x)  if (!(id_var %in% names(x)))
        stop("'id_var' must be included as a grouping factor in each ",
             "of the longitudinal submodels", call. = FALSE)) 
      return(id_var)
    }      
  } else {  # only one grouping factor (assumed to be subject ID)
    only_cnm <- unique(sapply(y_cnms, names))
    if (length(only_cnm) > 1L)
      stop("The grouping factor (ie, subject ID variable) is not the ",
           "same in all longitudinal submodels", call. = FALSE)
    if ((!is.null(id_var)) && (!identical(id_var, only_cnm)))
      warning("The user specified 'id_var' (", paste(id_var), 
              ") and the assumed ID variable based on the single ",
              "grouping factor (", paste(only_cnm), ") are not the same; ", 
              "'id_var' will be ignored", call. = FALSE, immediate. = TRUE)
    return(only_cnm)
  }
}


# Check the factor list corresponding to subject ID is the same in each 
# of the longitudinal submodels
#
check_id_list <- function(id_var, y_flist) {
  id_list <- unique(lapply(y_flist, function(x) unique(x[[id_var]])))
  if (length(id_list) > 1L)
    stop("The subject IDs are not the same in all longitudinal submodels.",
         call. = FALSE)
  unlist(id_list)  
}


# Function to check that the assoc argument only includes supported association
# types. The function returns a list with logicals specifying which association
# type have been requested.
# 
# @param m Integer specifying the longitudinal submodel 
# @param x The assoc argument specified by the user -- should be a character 
#   vector or NULL
# @param y_cnms The RE component names for each of the longitudinal submodels 
# @param supported_assocs A character vector showing the supported
#   association types
# @param id_var The name of the id variable 
# @return A list of logicals indicating the desired association types
validate_assoc <- function(m, x, y_cnms, supported_assocs, id_var, xmat) {
  
  # Select association structure for submodel m only
  # 'x_tmp' will be changed in code below, whilst 'x' stays equal to user input
  x_tmp <- x <- x[[m]]
  
  # Identify which association types were specified
  if (!is.null(x_tmp)) {
    x_tmp <- gsub("^etalag\\(.*", "etalag", x_tmp) 
    x_tmp <- gsub("^mulag\\(.*", "mulag", x_tmp) 
    x_tmp <- gsub("^shared_b.*", "shared_b", x_tmp) 
    x_tmp <- gsub("^shared_coef\\(.*", "shared_coef", x_tmp) 
    x_tmp <- gsub("^etavalue_interact\\(.*", "etavalue_interact", x_tmp) 
    x_tmp <- gsub("^muvalue_interact\\(.*", "muvalue_interact", x_tmp) 
    x_tmp <- gsub("^etaslope_interact\\(.*", "etaslope_interact", x_tmp) 
    x_tmp <- gsub("^muslope_interact\\(.*", "muslope_interact", x_tmp) 
  }
  assoc <- sapply(supported_assocs, function(y) y %in% x_tmp, simplify = FALSE)
  if (is.null(x_tmp)) {
    assoc$null <- TRUE
  } else if (is.character(x_tmp)) {
    if (!all(x_tmp %in% supported_assocs))
      stop("An unsupported association type has been specified. The ",
           "'assoc' argument can only include the following association ", 
           "types: ", paste(supported_assocs, collapse = ", "), call. = FALSE)
    if ((assoc$null) && (length(x_tmp) > 1L))
      stop("In 'assoc' argument, 'null' cannot be specified in ",
           "conjuction with another association type", call. = FALSE)
    if (assoc$etavalue && assoc$muvalue)
      stop("In 'assoc' argument, 'etavalue' and 'muvalue' cannot be specified ",
           "together", call. = FALSE)
    if (assoc$etaslope && assoc$muslope)
      stop("In 'assoc' argument, 'etaslope' and 'muslope' cannot be specified ",
           "together", call. = FALSE)   
    if (assoc$etavalue && assoc$muvalue)
      stop("In 'assoc' argument, 'etalag' and 'mulag' cannot be specified ",
           "together", call. = FALSE)
  } else { 
    stop("'assoc' argument should be a character vector or, for a multivariate ",
         "joint model, possibly a list of character vectors.", call. = FALSE)
  }
  
  # Identify which lags were requested
  if (assoc$etalag || assoc$mulag) {
    if (assoc$etalag) {
      val_etalag <- grep("^etalag.*", x, value = TRUE)
      val_lag <- unlist(strsplit(val_etalag, "etalag"))[-1]
    } else if (assoc$mulag) {
      val_mulag <- grep("^mulag.*", x, value = TRUE)
      val_lag <- unlist(strsplit(val_mulag, "mulag"))[-1]
    }
    if (length(val_lag)) {
      assoc$which_lag <- tryCatch(eval(parse(text = paste0("c", val_lag))), 
                                  error = function(x) 
                                    stop("Incorrect specification of the lagged ",
                                         "association structure. See Examples in help ",
                                         "file.", call. = FALSE))
      if (length(assoc$which_lag) > 1L) 
        stop("Currently only one lag time is allowed for the lagged association ",
             "structure.", call. = FALSE)
    } else {
      stop("'etalag' association structure was specified incorrectly. It should ",
           "include a suffix with the desired lag inside parentheses. See the ",
           "help file for details.", call. = FALSE)    
    }    
  } else assoc$which_lag <- 0  
  
  # Check interaction formula was specified correctly
  if (any(assoc$etavalue_interact, assoc$muvalue_interact, 
          assoc$etaslope_interact, assoc$muslope_interact)) {
    assoc$formula_interact <- sapply(c("etavalue_interact", "muvalue_interact",
                                       "etaslope_interact", "muslope_interact"),
                                     function(y, x) {
                                       if (assoc[[y]]) {
                                         val <- grep(paste0("^", y, ".*"), x, value = TRUE)
                                         val <- unlist(strsplit(val, y))[-1]
                                         fm <- tryCatch(eval(parse(text = val)), 
                                                        error = function(e) 
                                                          stop(paste0("Incorrect specification of the formula in the '", y,
                                                                      "' association structure. See Examples in help ",
                                                                      "file."), call. = FALSE))
                                         if (!is(fm, "formula"))
                                           stop(paste0("Suffix to '", y, "' association structure should include ",
                                                       "a formula within parentheses."), call. = FALSE)
                                         if (identical(length(fm), 3L))
                                           stop(paste0("Formula specified for '", y, "' association structure should not ",
                                                       "include a response."), call. = FALSE)
                                         if (length(lme4::findbars(fm)))
                                           stop(paste0("Formula specified for '", y, "' association structure should only ",
                                                       "include fixed effects."), call. = FALSE)
                                         if (fm[[2L]] == 1)
                                           stop(paste0("Formula specified for '", y, "' association structure cannot ",
                                                       "be an intercept only."), call. = FALSE)
                                         fm
                                       } else NULL
                                     }, x = x, simplify = FALSE, USE.NAMES = TRUE)
  } else assoc$formula_interact <- NULL  # no interaction assoc at all for submodel m
  
  # Identify which subset of shared random effects were specified
  max_which_b_zindex <- length(y_cnms[[m]][[id_var]])
  val_b <- grep("^shared_b.*", x, value = TRUE)
  val_coef <- grep("^shared_coef.*", x, value = TRUE)
  
  if (length(val_b)) {
    val_b <- unlist(strsplit(val_b, "shared_b"))[-1]
    if (length(val_b)) {
      assoc$which_b_zindex <- tryCatch(eval(parse(text = paste0("c", val_b))), 
                                       error = function(x) 
                                         stop("Incorrect specification of the 'shared_b' ",
                                              "association structure. See Examples in help ",
                                              "file.", call. = FALSE))
    } else assoc$which_b_zindex <- seq_len(max_which_b_zindex)
    if (any(assoc$which_b_zindex > max_which_b_zindex))
      stop(paste0("The indices specified for the shared random effects association ",
                  "structure are greater than the number of subject-specific random ", 
                  "effects (this error was encountered for longitudinal submodel ", m,
                  ")"), call. = FALSE)
    names(assoc$which_b_zindex) <- y_cnms[[m]][[id_var]][assoc$which_b_zindex]
  } else assoc$which_b_zindex <- numeric(0)
  
  if (length(val_coef)) {
    val_coef <- unlist(strsplit(val_coef, "shared_coef"))[-1]
    if (length(val_coef)) {
      assoc$which_coef_zindex <- tryCatch(eval(parse(text = paste0("c", val_coef))), 
                                          error = function(x) 
                                            stop("Incorrect specification of the 'shared_coef' ",
                                                 "association structure. See Examples in help ",
                                                 "file.", call. = FALSE))
    } else assoc$which_coef_zindex <- seq_len(max_which_b_zindex)
    if (any(assoc$which_coef_zindex > max_which_b_zindex))
      stop(paste0("The indices specified for the shared random effects association ",
                  "structure are greater than the number of subject-specific random ", 
                  "effects (this error was encountered for longitudinal submodel ", m,
                  ")"), call. = FALSE)
  } else assoc$which_coef_zindex <- numeric(0)  
  
  if (length(intersect(assoc$which_b_zindex, assoc$which_coef_zindex)))
    stop("The same random effects indices should not be specified in both ",
         "'shared_b' and 'shared_coef'. Specifying indices in 'shared_coef' ",
         "will include both the fixed and random components.", call. = FALSE)
  
  if (length(assoc$which_coef_zindex)) {
    if (length(y_cnms[[m]]) > 1L)
      stop("'shared_coef' association structure cannot be used when there is ",
           "clustering at levels other than the individual-level.", call. = FALSE)
    b_nms <- y_cnms[[m]][[id_var]][assoc$which_coef_zindex]
    names(assoc$which_coef_zindex) <- b_nms
    beta_nms <- colnames(xmat[[m]])
    assoc$which_coef_xindex <- sapply(b_nms, function(x, beta_nms) {
      beta_match <- grep(x, beta_nms, fixed = TRUE)
      if (!length(beta_match)) {
        stop("In association structure 'shared_coef', no matching fixed effect ",
             "component could be found for the following random effect: ", x, 
             ". Perhaps consider using 'shared_b' association structure instead.")
      } else if (length(beta_match) > 1L) {
        stop("Bug found: In association structure 'shared_coef', multiple ",
             "fixed effect components have been found to match the following ",
             "random effect: ", x)
      }  
      beta_match
    }, beta_nms)
  } else {
    assoc$which_coef_xindex <- numeric(0)
  }
  
  if (!identical(length(assoc$which_coef_zindex), length(assoc$which_coef_xindex)))
    stop("Bug found: the lengths of the fixed and random components of the ",
         "'shared_coef' association structure are not the same.")

  assoc
}


# Function to return the variables used in fitting a model, as well as the
# predvars, for both the fixed and random parts of a (g)lmer fit
#
# @param mod A (g)lmer model object
# @return A named list of lists
get_formvars <- function(mod) {
  vars_f <- grep("", attr(terms(mod, fixed.only = TRUE), "variables"), value = TRUE)
  vars_r <- grep("", attr(terms(mod, random.only = TRUE), "variables"), value = TRUE)
  predvars_f <- grep("", attr(terms(mod, fixed.only = TRUE), "predvars"), value = TRUE)
  predvars_r <- grep("", attr(terms(mod, random.only = TRUE), "predvars"), value = TRUE)
  list(formvars = list(fixed = vars_f, random = vars_r),
       predvars = list(fixed = predvars_f, random = predvars_r))
}


# Function to substitute variables in the formula of a fitted model
#
# @param mod A (g)lmer model object from which to extract the model formula
# @param formvars A list of the original variables used in the model formula
# @param newvars A list of the variables to use as replacements
# @return The reformulated model formula with the variables in oldvars replaced 
#   the corresponding entries in newvars
use_these_vars <- function(mod, oldvars, newvars) {
  if (!identical(length(oldvars), length(newvars)))
    stop("oldvars and newvars should be the same length.")
  fm <- formula(mod)
  if (!identical(oldvars, newvars)) {
    for (j in 1:length(oldvars))
      fm <- gsub(oldvars[[j]], newvars[[j]], fm, fixed = TRUE)    
    fm <- reformulate(fm[[3]], response = fm[[2]])
  }
  fm
}

# Function to substitute variables in the formula of a fitted model
# with the corresponding predvars based on the terms object for the model
#
# @param formula A formula object from a fitted model
# @param terms A terms object from a fitted model
# @return A reformulated model formula with variables replaced by predvars
use_predvars <- function(mod) {
  fm <- formula(mod)
  ff <- grep("", attr(terms(mod, fixed.only = TRUE), "variables"), value = TRUE)[-1]
  fr <- grep("", attr(terms(mod, random.only = TRUE), "variables"), value = TRUE)[-1]
  pf <- grep("", attr(terms(mod, fixed.only = TRUE), "predvars"), value = TRUE)[-1]
  pr <- grep("", attr(terms(mod, random.only = TRUE), "predvars"), value = TRUE)[-1]
  if (!identical(c(ff, fr), c(pf, pr))) {
    for (j in 1:length(ff))
      fm <- gsub(ff[[j]], pf[[j]], fm, fixed = TRUE)    
    for (j in 1:length(fr))
      fm <- gsub(fr[[j]], pr[[j]], fm, fixed = TRUE)    
    fm <- reformulate(fm[[3]], response = fm[[2]])
  }
  fm
}


# Function to return control arguments for lmer or glmer call
#
get_control_args <- function(glmer = FALSE, norank = FALSE) {
  if (glmer) {
    if (norank) {
      lme4::glmerControl(check.nlev.gtreq.5 = "ignore",
                         check.nlev.gtr.1 = "stop",
                         check.nobs.vs.rankZ = "ignore",
                         check.nobs.vs.nlev = "ignore",
                         check.nobs.vs.nRE = "ignore",
                         check.rankX = "ignore")      
    } else {
      lme4::glmerControl(check.nlev.gtreq.5 = "ignore",
                         check.nlev.gtr.1 = "stop",
                         check.nobs.vs.rankZ = "ignore",
                         check.nobs.vs.nlev = "ignore",
                         check.nobs.vs.nRE = "ignore")      
    }    
  } else {
    if (norank) {
      lme4::lmerControl(check.nlev.gtreq.5 = "ignore",
                        check.nlev.gtr.1 = "stop",
                        check.nobs.vs.rankZ = "ignore",
                        check.nobs.vs.nlev = "ignore",
                        check.nobs.vs.nRE = "ignore",
                        check.rankX = "ignore")       
    } else {
      lme4::lmerControl(check.nlev.gtreq.5 = "ignore",
                        check.nlev.gtr.1 = "stop",
                        check.nobs.vs.rankZ = "ignore",
                        check.nobs.vs.nlev = "ignore",
                        check.nobs.vs.nRE = "ignore")      
    }
  }
}

# Function to return standardised GK quadrature points and weights
#
# @return A list with two named elements (points and weights) each
#   of which is a numeric vector with length equal to the number of
#   quadrature nodes
get_quadpoints <- function(nodes = 15) {
  if (!is.numeric(nodes) || (length(nodes) > 1L)) {
    stop("'quadnodes' should be a numeric vector of length 1.")
  } else if (nodes == 15) {
    list(
      points = c(
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
        0.991455371120812639207),
      weights = c(
        0.0229353220105292249637,
        0.063092092629978553291,
        0.10479001032225018384,
        0.140653259715525918745,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.204432940075298892414,
        0.209482141084727828013,
        0.204432940075298892414,
        0.1903505780647854099133,
        0.169004726639267902827,
        0.140653259715525918745,
        0.1047900103222501838399,
        0.063092092629978553291,
        0.0229353220105292249637))      
  } else if (nodes == 11) {
    list(
      points = c(
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
        0.984085360094842464496),
      weights = c(
        0.042582036751081832865,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.241040339228647586701,
        0.186800796556492657467,
        0.115233316622473394025,
        0.042582036751081832865))     
  } else if (nodes == 7) {
    list(
      points = c(
        -0.9604912687080202834235,
        -0.7745966692414833770359,
        -0.4342437493468025580021,
        0,
        0.4342437493468025580021,
        0.7745966692414833770359,
        0.9604912687080202834235),
      weights = c(
        0.1046562260264672651938,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194))      
  } else stop("'quadnodes' must be either 7, 11 or 15.")  
}


# Function to return the appropriate intercept type based on family and link
#
# @param family A GLM family
# @param link A link function
# @return A list specifying whether an unbounded, lower bounded, or upper 
#   bounded intercept is required
check_intercept <- function(family, link) {
  if (family == "binomial") {
    if (link == "log") {
      unbound <- 0L; lobound <- 0L; upbound <- 1L;  # binomial, log
    } else {
      unbound <- 1L; lobound <- 0L; upbound <- 0L;  # binomial, !log
    }
  } else {
    if (link == "log") {
      unbound <- 1L; lobound <- 0L; upbound <- 0L;  # gamma/inv-gaus/poisson/nb, log
    } else {
      if (family == "gaussian") {
        unbound <- 1L; lobound <- 0L; upbound <- 0L;  # gaussian, any link
      } else {
        unbound <- 0L; lobound <- 1L; upbound <- 0L;  # gamma/inv-gaus/poisson/nb, !log      
      }
    }
  }  
  nlist(unbound, lobound, upbound) 
}


# Function to check if the response vector is real or integer
#
# @param family A GLM family
# @return A logical specify whether the response is real (TRUE) or integer (FALSE)
check_response_real <- function(family) {
  if (family == "binomial" || 
      family == "poisson" ||
      family == "neg_binomial_2") {
    FALSE
  } else TRUE
}


# Function to check if the submodel should include a dispersion term
#
# @param family A GLM family
# @return A logical specify whether the submodel includes a dispersion term
check_for_dispersion <- function(family) {
  if (family == "binomial" || family == "poisson") FALSE else TRUE
}


# Function to calculate the number of association parameters in the model
#
# @param has_assoc A named list specifying whether each longitudinal submodel 
#   is linked to the event outcome using each potential type of association structure
# @param which_b_zindex A list of numeric vectors indicating the random effects from each
#   longitudinal submodel that are to be used in the shared_b association structure
# @param which_coef_zindex A list of numeric vectors indicating the random effects from each
#   longitudinal submodel that are to be used in the shared_coef association structure
# @return Integer indicating the number of association parameters in the model 
get_num_assoc_pars <- function(has_assoc, which_b_zindex, which_coef_zindex) {
  sel <- c("etavalue", "etaslope", "etalag", "etaauc", "muvalue", "muslope", "mulag", "muauc")
  a_K <- sum(unlist(has_assoc[sel]))
  a_K <- a_K + length(unlist(which_b_zindex)) + length(unlist(which_coef_zindex))
  return(a_K)
}


#-------- The following functions could be moved to misc.R

# Remove a specified character string from the names of an
# object (for example, a matched call)
#
# @param x The matched call
# @param string The character string to be removed
strip_nms <- function(x, string) {
  names(x) <- gsub(string, "", names(x))
  x
}

# Supplies names for the output list returned by most stanjm methods
#
# @param x The list object to which the names are to be applied
# @param M The number of longitudinal submodels
list_nms <- function(object, M) {
  if (!is.list(object)) stop("'object' argument should be a list")
  nms <- paste0("Long", 1:M)
  if (length(object) > M) nms <- c(nms, "Event")
  names(object) <- nms
  object
}


#-------- The following functions could be moved to priors.R

#' @rdname priors
#' @export 
#' @param prior_scale_for_dispersion Prior scale for the standard error of the 
#'   regression in Gaussian models, which is given a half-Cauchy prior truncated
#'   at zero.
#' @param min_prior_scale Minimum prior scale for the intercept and 
#'   coefficients.
#' @param scaled A logical scalar, defaulting to \code{TRUE}. If \code{TRUE} the
#'   \code{prior_scale} is further scaled by the range of the predictor if the 
#'   predictor has exactly two unique values and scaled by twice the standard
#'   deviation of the predictor if it has more than two unique values.
#'
priorLong_options <- function(prior_scale_for_dispersion = 5, 
                              min_prior_scale = 1e-12, 
                              scaled = TRUE) {
  validate_parameter_value(prior_scale_for_dispersion)
  validate_parameter_value(min_prior_scale)
  nlist(scaled, min_prior_scale, prior_scale_for_dispersion)
}

#' @rdname priors
#' @export 
#' @param prior_scale_for_basehaz Usage depends on the baseline hazard 
#'   specified in the \code{base_haz} argument of the  
#'   \code{\link{stan_jm}} call. 
#'   If \code{base_haz = "weibull"} then this argument specifies the
#'   prior scale for the shape parameter of the Weibull distribution, 
#'   which is given a half-Cauchy truncated at zero.
#'   If \code{base_haz = "splines"} then this argument specifies the
#'   prior scale(s) for the spline regression coefficients, 
#'   which are given normal distributions with location 0.
#'   If \code{base_haz = "piecewise"} then this argument specifies the
#'   prior scale for the parameters corresponding to the log baseline 
#'   hazard within each interval, which are given normal distributions
#'   with location 0.
#'  
priorEvent_options <- function(prior_scale_for_basehaz = 50,
                               min_prior_scale = 1e-12, 
                               scaled = TRUE) {
  validate_parameter_value(prior_scale_for_basehaz)
  validate_parameter_value(min_prior_scale)
  nlist(scaled, min_prior_scale, prior_scale_for_basehaz)
}

#' @rdname priors
#' @export 
#'
priorAssoc_options <- function(min_prior_scale = 1e-12) {
  validate_parameter_value(min_prior_scale)
  nlist(min_prior_scale)
}



