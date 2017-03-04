# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Monash University
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
#' @export
#' @templateVar rareargs na.action,contrasts
#' @template args-same-as-rarely
#' @template args-dots
#' @template args-prior_covariance
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template args-sparse
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
#' @param offset Not currently implemented. Same as \code{\link[stats]{glm}}.
#' @param family The family (and possibly also the link function) for the 
#'   longitudinal submodel(s). See \code{\link[lme4]{glmer}} for details. 
#'   If fitting a multivariate joint model, then this can optionally be a
#'   list of families, in which case each element of the list specifies the
#'   family for one of the longitudinal submodels.
#' @param assoc A character string or character vector specifying the joint
#'   model association structure. Possible association structures that can
#'   be used include: "etavalue" (the default); "etaslope"; "etalag"; "etaauc"; 
#'   "muvalue"; "muslope"; "mulag"; "muauc"; "shared_b"; "shared_coef"; or "null". 
#'   These are described in the \strong{Details} section below. For a multivariate 
#'   joint model, different association structures can optionally be used for 
#'   each longitudinal submodel by specifying a list of character
#'   vectors, with each element of the list specifying the desired association 
#'   structure for one of the longitudinal submodels. Specifying \code{assoc = NULL}
#'   will fit a joint model with no association structure (equivalent  
#'   to fitting separate longitudinal and time-to-event models). It is also 
#'   possible to include interaction terms between the association term 
#'   ("etavalue", "etaslope", "etalag", "etaauc, "muvalue", "muslope", "mulag",
#'   or "muauc") with observed data/covariates. It is also possible, when fitting 
#'   a multivariate joint model, to include interaction terms between the 
#'   association terms ("etavalue" or "muvalue") corresponding to the different
#'   longitudinal outcomes. See the \strong{Details} section as well as the
#'   \strong{Examples} below.
#' @param basehaz A character string indicating which baseline hazard to use
#'   for the event submodel. Options are a Weibull baseline hazard
#'   (\code{"weibull"}, the default), a B-splines approximation estimated 
#'   for the log baseline hazard (\code{"bs"}), or a piecewise
#'   constant baseline hazard (\code{"piecewise"}).
#' @param basehaz_ops A named list specifying options related to the baseline
#'   hazard. Currently this can include: \cr
#'   \describe{
#'     \item{\code{df}}{A positive integer specifying the degrees of freedom 
#'     for the B-splines if \code{basehaz = "bs"}, or the number of
#'     intervals used for the piecewise constant baseline hazard if 
#'     \code{basehaz = "piecewise"}. The default is 6.}
#'     \item{\code{knots}}{An optional numeric vector specifying the internal knot 
#'     locations for the B-splines if \code{basehaz = "bs"}, or the 
#'     internal cut-points for defining intervals of the piecewise constant 
#'     baseline hazard if \code{basehaz = "piecewise"}. Knots cannot be
#'     specified if \code{df} is specified. If not specified, then the 
#'     default is to use \code{df - 4} knots if \code{basehaz = "bs"},
#'     or \code{df - 1} knots if \code{basehaz = "piecewise"}, which are
#'     placed at equally spaced percentiles of the distribution of
#'     observed event times.}
#'   }
#' @param dataAssoc A data frame containing observed covariates used in the
#'   interactions between association terms and observed data. See the
#'   \strong{Details} and \strong{Examples} sections for details on how to 
#'   specify the formulas for the interactions as part of the \code{assoc}  
#'   argument. If interactions between association terms and observed data are
#'   specified as part of the \code{assoc} argument, but \code{dataAssoc} is 
#'   not specified, then the values for the covariates will be taken from the 
#'   data frame(s) provided in \code{dataLong}. (Specifying \code{dataAssoc}
#'   directly means that a different measurement time schedule could be used 
#'   for the covariates in the association term interactions when compared with
#'   the measurement time schedule used for the longitudinal outcomes and  
#'   covariates provided in \code{dataLong}.
#' @param quadnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard in the likelihood function. 
#'   Options are 15 (the default), 11 or 7.
#' @param subsetLong,subsetEvent Same as subset in \code{\link[stats]{glm}}.
#'   However, if fitting a multivariate joint model and a list of data frames 
#'   is provided in \code{dataLong} then a corresponding list of subsets 
#'   must be provided in \code{subsetLong}.
#' @param weights Experimental and should be used with caution. The 
#'   user can optionally supply a 2-column data frame containing a set of
#'   'prior weights' to be used in the estimation process. The data frame should
#'   contain two columns: the first containing the IDs for each individual, and 
#'   the second containing the corresponding weights. The data frame should only
#'   have one row for each individual; that is, weights should be constant 
#'   within individuals.
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
#' @param priorLong,priorEvent,priorAssoc The prior distributions for the 
#'   regression coefficients in the longitudinal submodel(s), event submodel,
#'   and the association parameter(s).
#'    
#'   Can be a call to one of the various functions provided by 
#'   \pkg{rstanarm} for specifying priors. The subset of these functions that 
#'   can be used for the prior on the coefficients can be grouped into several 
#'   "families":
#'   
#'   \tabular{ll}{
#'     \strong{Family} \tab \strong{Functions} \cr 
#'     \emph{Student t family} \tab \code{normal}, \code{student_t}, \code{cauchy} \cr 
#'     \emph{Hierarchical shrinkage family} \tab \code{hs}, \code{hs_plus} \cr 
#'     \emph{Laplace family} \tab \code{laplace}, \code{lasso} \cr
#'     \emph{Product normal family} \tab \code{product_normal} \cr
#'   }
#'   
#'   See the \link[=priors]{priors help page} for details on the families and 
#'   how to specify the arguments for all of the functions in the table above.
#'   To omit a prior ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior} can be set to \code{NULL}, although this is rarely a good
#'   idea.
#'   
#'   \strong{Note:} Unless \code{QR=TRUE}, if \code{prior} is from the Student t
#'   family or Laplace family, and if the \code{autoscale} argument to the 
#'   function used to specify the prior (e.g. \code{\link{normal}}) is left at 
#'   its default and recommended value of \code{TRUE}, then the default or 
#'   user-specified prior scale(s) may be adjusted internally based on the scales
#'   of the predictors. See the \link[=priors]{priors help page} for details on
#'   the rescaling and the \code{\link{prior_summary}} function for a summary of
#'   the priors used for a particular model.
#' @param priorLong_intercept,priorEvent_intercept The prior distributions  
#'   for the intercepts in the longitudinal submodel(s) and event submodel. 
#'   Can be a call to \code{normal}, \code{student_t} or 
#'   \code{cauchy}. See the \link[=priors]{priors help page} for details on 
#'   these functions. To omit a prior on the intercept ---i.e., to use a flat
#'   (improper) uniform prior--- \code{prior_intercept} can be set to
#'   \code{NULL}.
#'   
#'   \strong{Note:} If using a dense representation of the design matrix 
#'   ---i.e., if the \code{sparse} argument is left at its default value of
#'   \code{FALSE}--- then the prior distribution for the intercept is set so it
#'   applies to the value when all predictors are centered.
#' @param priorLong_aux The prior distribution for the "auxiliary" parameters
#'   in the longitudinal submodels (if applicable). 
#'   The "auxiliary" parameter refers to a different parameter 
#'   depending on the \code{family}. For Gaussian models \code{priorLong_aux} 
#'   controls \code{"sigma"}, the error 
#'   standard deviation. For negative binomial models \code{priorLong_aux} controls 
#'   \code{"reciprocal_dispersion"}, which is similar to the 
#'   \code{"size"} parameter of \code{\link[stats]{rnbinom}}:
#'   smaller values of \code{"reciprocal_dispersion"} correspond to 
#'   greater dispersion. For gamma models \code{priorLong_aux} sets the prior on 
#'   to the \code{"shape"} parameter (see e.g., 
#'   \code{\link[stats]{rgamma}}), and for inverse-Gaussian models it is the 
#'   so-called \code{"lambda"} parameter (which is essentially the reciprocal of
#'   a scale parameter). Binomial and Poisson models do not have auxiliary 
#'   parameters. 
#'   
#'   \code{priorLong_aux} can be a call to \code{exponential} to 
#'   use an exponential distribution, or \code{normal}, \code{student_t} or 
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy 
#'   prior. See \code{\link{priors}} for details on these functions. To omit a 
#'   prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{priorLong_aux} to \code{NULL}.
#'   
#'   If fitting a multivariate joint model, you have the option to
#'   specify a list of prior distributions, however the elements of the list
#'   that correspond to any longitudinal submodel which does not have an 
#'   auxiliary parameter will be ignored. 
#' @param priorEvent_aux The prior distribution for the "auxiliary" parameters
#'   in the event submodel. The "auxiliary" parameters refers to different  
#'   parameters depending on the baseline hazard. For \code{basehaz = "weibull"}
#'   the auxiliary parameter is the Weibull shape parameter. For 
#'   \code{basehaz = "bs"} the auxiliary parameters are the coefficients for the
#'   B-spline approximation to the log baseline hazard.
#'   For \code{basehaz = "piecewise"} the auxiliary parameters are the piecewise
#'   estimates of the log baseline hazard. 
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
#'   using either a Weibull distribution (\code{basehaz = "weibull"}) or a
#'   piecewise constant baseline hazard (\code{basehaz = "piecewise"}), or 
#'   approximated using cubic B-splines (\code{basehaz = "bs"}). 
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
#'               assoc = c("etavalue", "etavalue_data(~ trt)"))
#' summary(f4)  
#'                     
#' }
#' 
#' @import data.table
#' @importFrom lme4 lmerControl glmerControl
#' 
stan_jm <- function(formulaLong, dataLong, formulaEvent, dataEvent, time_var, 
                    id_var, family = gaussian, assoc = "etavalue", dataAssoc,
                    basehaz = c("weibull", "bs", "piecewise"), basehaz_ops, 
                    quadnodes = 15, subsetLong, subsetEvent, init = "model_based", 
                    na.action = getOption("na.action", "na.omit"), weights, 
                    offset, contrasts, ...,				          
                    priorLong = normal(), priorLong_intercept = normal(), 
                    priorLong_aux = cauchy(0, 5), priorEvent = normal(), 
                    priorEvent_intercept = normal(), priorEvent_aux = cauchy(0, 50),
                    priorAssoc = normal(), prior_covariance = decov(), prior_PD = FALSE, 
                    algorithm = c("sampling", "meanfield", "fullrank"), 
                    adapt_delta = NULL, max_treedepth = NULL, QR = FALSE, 
                    sparse = FALSE, long_lp = TRUE, event_lp = TRUE) {
  
  
  #=============================
  # Pre-processing of arguments
  #=============================  
  
  # Check for arguments not yet implemented
  if (!missing(offset)) 
    stop("Offsets are not yet implemented for stan_jm")
  if (QR)               
    stop("QR decomposition not yet implemented for stan_jm")
  if (sparse)
    stop("'sparse' option is not yet implemented for stan_jm")
  #  if (algorithm %in% c("meanfield", "fullrank"))
  #    stop ("Meanfield and fullrank algorithms not yet implemented for stan_jm")
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  if (missing(weights))     weights     <- NULL
  if (missing(id_var))      id_var      <- NULL
  if (missing(subsetLong))  subsetLong  <- NULL
  if (missing(dataAssoc))   dataAssoc   <- NULL
  
  basehaz   <- match.arg(basehaz)
  algorithm <- match.arg(algorithm)
  
  # Validate arguments
  formulaLong <- validate_arg(formulaLong, "formula")
  M           <- length(formulaLong)
  dataLong    <- validate_arg(dataLong,   "data.frame", null_ok = TRUE, validate_length = M)
  subsetLong  <- validate_arg(subsetLong, "vector",     null_ok = TRUE, validate_length = M)
  assoc       <- validate_arg(assoc,      "character",  null_ok = TRUE, validate_length = M, broadcast = TRUE)

  # Check family and link
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  if (!is.list(family)) {
    family <- rep(list(family), M) 
  } else if (!length(family) == M) {
    stop("family is a list of the incorrect length.")
  }
  family <- lapply(family, validate_family)
  fam <- lapply(family, function(x) 
    which(pmatch(supported_families, x$family, nomatch = 0L) == 1L))
  if (any(lapply(fam, length) == 0L)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  supported_links <- lapply(fam, function(x) supported_glm_links(supported_families[x]))
  link <- mapply(function(x, i) which(supported_links[[i]] == x$link),
                 family, seq_along(family), SIMPLIFY = TRUE)
  if (any(lapply(link, length) == 0L)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))

  # Matched call
  call <- match.call(expand.dots = TRUE)    
  mc   <- match.call(expand.dots = FALSE)
  mc$time_var <- mc$id_var <- mc$assoc <- 
    mc$basehaz <- mc$basehaz_ops <-
    mc$df <- mc$knots <- mc$quadnodes <- NULL
  mc$priorLong <- mc$priorLong_intercept <- mc$priorLong_aux <-
    mc$priorEvent <- mc$priorEvent_intercept <- mc$priorEvent_aux <-
    mc$priorAssoc <- mc$prior_covariance <-
    mc$prior_PD <- mc$algorithm <- mc$scale <- 
    mc$concentration <- mc$shape <- mc$init <- 
    mc$adapt_delta <- mc$max_treedepth <- 
    mc$... <- mc$QR <- NULL
  mc$weights <- NULL 
  mc$long_lp <- mc$event_lp <- NULL

  # Create call for longitudinal submodel  
  y_mc <- mc
  y_mc <- strip_nms(y_mc, "Long") 
  y_mc$formulaEvent <- y_mc$dataEvent <- y_mc$subsetEvent <- NULL

  # Create call for each longitudinal submodel separately
  m_mc <- lapply(1:M, function(m, old_call) {
    new_call <- old_call
    fm       <- eval(old_call$formula)
    data     <- eval(old_call$data)
    subset   <- eval(old_call$subset)
    family   <- eval(old_call$family)
    new_call$formula <- if (is(fm, "list"))     fm[[m]]     else old_call$formula
    new_call$data    <- if (is(data, "list"))   data[[m]]   else old_call$data
    new_call$subset  <- if (is(subset, "list")) subset[[m]] else old_call$subset
    new_call$family  <- if (is(family, "list")) family[[m]] else old_call$family
    new_call
  }, old_call = y_mc)

  # Create call for event submodel
  e_mc <- mc
  e_mc <- strip_nms(e_mc, "Event")
  e_mc$formulaLong <- e_mc$dataLong <- e_mc$family <- e_mc$subsetLong <- NULL
  
  # Is priorLong* already a list?
  priorLong           <- maybe_broadcast_priorarg(priorLong)
  priorLong_intercept <- maybe_broadcast_priorarg(priorLong_intercept)
  priorLong_aux       <- maybe_broadcast_priorarg(priorLong_aux)
    
  #================================
  # Data for longitudinal submodel
  #================================
  
  # Fit separate longitudinal submodels
  y_mod_stuff <- mapply(handle_glmod, m_mc, family, 
                        MoreArgs = list(supported_families = supported_families, 
                                        supported_links = supported_links, 
                                        sparse = sparse), SIMPLIFY = FALSE)

  # Construct single cnms list for all longitudinal submodels
  cnms <- get_common_cnms(fetch(y_mod_stuff, "cnms"))
  cnms_nms <- names(cnms)
  
  # Additional error checks
  id_var <- check_id_var (id_var, fetch(y_mod_stuff, "cnms"))
  unique_id_list <- check_id_list(id_var, fetch(y_mod_stuff, "flist"))
  
  # Construct prior weights
  has_weights <- (!is.null(weights))
  if (has_weights) check_arg_weights(weights, id_var)
  y_weights <- lapply(y_mod_stuff, handle_weights, weights, id_var)
  
  #=========================
  # Data for event submodel
  #=========================
  
  if (!id_var %in% colnames(dataEvent))
    stop(paste0("Variable '", id_var, "' must be appear in dataEvent"), call. = FALSE)
  
  # Fit separate event submodel
  e_mod_stuff <- handle_coxmod(e_mc, quadnodes = quadnodes, id_var = id_var, 
                                unique_id_list = unique_id_list, sparse = sparse)
  
  # Construct prior weights
  e_weights <- handle_weights(e_mod_stuff, weights, id_var)

  # Baseline hazard
  ok_basehaz <- nlist("weibull", "bs", "piecewise")
  basehaz <- handle_basehaz(basehaz, basehaz_ops, ok_basehaz = ok_basehaz, 
                            eventtime = e_mod_stuff$eventtime, d = e_mod_stuff$d)
  
  # Incorporate intercept term if Weibull baseline hazard
  e_has_intercept <- e_mod_stuff$has_intercept <- (basehaz$type == 1L)

  #================================
  # Data for association structure
  #================================
  
  # Handle association structure
  ok_assoc <- c("null", "etavalue", "muvalue", "etaslope", "muslope", 
                "etalag", "mulag", "etaauc", "muauc", "shared_b", "shared_coef")
  ok_assoc_data         <- ok_assoc[2:5]
  ok_assoc_interactions <- ok_assoc[2:3]
  
  assoc <- mapply(validate_assoc, assoc, y_mod_stuff, 
                  MoreArgs = list(ok_assoc = ok_assoc, ok_assoc_data = ok_assoc_data,
                                  ok_assoc_interactions = ok_assoc_interactions, 
                                  id_var = id_var, M = M))
  assoc <- check_order_of_assoc_interactions(assoc, ok_assoc_interactions)
  colnames(assoc) <- paste0("Long", 1:M)
  
  # Indicator of each association type, for each longitudinal submodel
  sel <- grep("which_", rownames(assoc), invert = TRUE)
  has_assoc <- as.integer(assoc[sel,])  # Integer instead of logical and no which_{*} information

  # Time shift used for numerically calculating derivative of linear predictor 
  # or expected value of longitudinal outcome using one-sided difference
  eps <- 1E-5
  
  # Unstandardised quadrature nodes for AUC association structure
  auc_quadnodes <- 15
  auc_quadpoints <- get_quadpoints(auc_quadnodes)
  auc_quadweights <- unlist(
    lapply(e_mod_stuff$quadtimes, function(x) 
      lapply(x, function(y) 
        lapply(auc_quadpoints$weights, unstandardise_quadweights, 0, y))))

  # Return design matrices for evaluating longitudinal submodel quantities
  # at the quadrature points
  a_mod_stuff <- mapply(handle_a_mod, 1:M, m_mc, y_mod_stuff, SIMPLIFY = FALSE,
                        MoreArgs = list(e_mod_stuff = e_mod_stuff, assoc = assoc, 
                                        id_var = id_var, time_var = time_var, eps = eps))
  
  # Number of association parameters
  a_K <- get_num_assoc_pars(assoc, a_mod_stuff)

  #=====================
  # Prior distributions
  #=====================
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso", "product_normal")
  ok_intercept_dists <- ok_dists[1:3]
  ok_y_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  ok_e_aux_dists <- ok_dists[1:3]
  
  # Priors for longitudinal submodel(s)
  y_prior_stuff <- mapply(
    handle_glm_prior,
    priorLong,
    nvars = fetch(y_mod_stuff, "K"),
    link = fetch(family, "link"),
    MoreArgs = list(default_scale = 2.5, ok_dists = ok_dists), 
    SIMPLIFY = FALSE)

  if (length(unique(fetch(y_prior_stuff, "dist"))) > 1L)  # temporary
    stop("stan_jm does not currently allow different prior distributions to ",
         "be used for each of the longitudinal submodels.")
  
  y_prior_intercept_stuff <- mapply(
    handle_glm_prior,
    priorLong_intercept,
    link = fetch(family, "link"),
    MoreArgs = list(nvars = 1, default_scale = 10, ok_dists = ok_intercept_dists), 
    SIMPLIFY = FALSE)

  y_prior_aux_stuff <- mapply(
    handle_glm_prior,
    priorLong_aux,
    MoreArgs = list(nvars = 1, default_scale = 5, link = NULL, ok_dists = ok_y_aux_dists), 
    SIMPLIFY = FALSE)  
  
  # Priors for event submodel
  e_prior_stuff <- 
    handle_glm_prior(
      priorEvent,
      nvars = e_mod_stuff$K,
      default_scale = 2.5,
      link = NULL,
      ok_dists = ok_dists
    )  
  
  e_prior_intercept_stuff <- 
    handle_glm_prior(
      priorEvent_intercept,
      nvars = 1,
      default_scale = 50,
      link = NULL,
      ok_dists = ok_intercept_dists
    )  
  
  e_prior_aux_stuff <-
    handle_glm_prior(
      priorEvent_aux,
      nvars = basehaz$df,
      default_scale = 50,
      link = NULL,
      ok_dists = ok_e_aux_dists
    )
  
  # Priors for association parameters
  a_prior_stuff <- 
    handle_glm_prior(
      priorAssoc,
      nvars = a_K,
      default_scale = 25,
      link = NULL,
      ok_dists = ok_dists
    )  

  # Minimum scaling of priors
  y_prior_stuff           <- Map(autoscale_prior, y_prior_stuff, y_mod_stuff, QR = QR, use_x = TRUE)
  y_prior_intercept_stuff <- Map(autoscale_prior, y_prior_intercept_stuff, y_mod_stuff, QR = QR, use_x = FALSE)
  y_prior_aux_stuff       <- Map(autoscale_prior, y_prior_aux_stuff, y_mod_stuff, QR = QR, use_x = FALSE)
  e_prior_stuff           <- autoscale_prior(e_prior_stuff, e_mod_stuff, QR = QR, use_x = TRUE)
  e_prior_intercept_stuff <- autoscale_prior(e_prior_intercept_stuff, e_mod_stuff, QR = QR, use_x = FALSE)
  e_prior_aux_stuff       <- autoscale_prior(e_prior_aux_stuff, e_mod_stuff, QR = QR, use_x = FALSE)
  a_prior_stuff           <- autoscale_prior(a_prior_stuff, e_mod_stuff, QR = QR, use_x = FALSE)
    
  #=========================
  # Data for export to Stan
  #=========================
  
  standata <- list(  
    # dimensions
    quadnodes    = as.integer(quadnodes),
    has_weights  = as.integer(has_weights),
    M            = as.integer(M),
    Npat         = as.integer(e_mod_stuff$Npat),
    e_K          = as.integer(e_mod_stuff$K),
    a_K          = as.integer(a_K),
    y_K          = fetch_array(y_mod_stuff, "K"), 
    y_N          = fetch_array(y_mod_stuff, "N"), 
    y_real_N     = fetch_array(y_mod_stuff, "real_N"), 
    y_int_N      = fetch_array(y_mod_stuff, "int_N"), 
    y_N01        = as.array(t(sapply(fetch(y_mod_stuff, "N01"), cbind))), 
    Npat_times_quadnodes = as.integer(e_mod_stuff$Npat * quadnodes),
    
    # data for longitudinal submodel(s)
    link = as.array(link),
    y_has_intercept         = as.array(as.integer(fetch_(y_mod_stuff, "has_intercept"))),
    y_has_intercept_unbound = as.array(as.integer(fetch_(y_mod_stuff, "has_intercept_unbound"))),
    y_has_intercept_lobound = as.array(as.integer(fetch_(y_mod_stuff, "has_intercept_lobound"))),
    y_has_intercept_upbound = as.array(as.integer(fetch_(y_mod_stuff, "has_intercept_upbound"))),
    y_has_aux               = as.array(as.integer(fetch_(y_mod_stuff, "has_aux"))),
    y_xbar                  = fetch_array(y_mod_stuff, "xbar"),
    trials                  = fetch_array(y_mod_stuff, "trials"),
    y_weights               = as.array(unlist(y_weights)),
    y_X = as.array(as.matrix(Matrix::bdiag(fetch(y_mod_stuff, "xtemp")))),
    
    # data for event submodel
    basehaz_type    = as.integer(basehaz$type),
    basehaz_df      = as.integer(basehaz$df),
    e_has_intercept = as.integer(e_has_intercept),
    nrow_y_Xq       = NROW(a_mod_stuff[[1]]$xtemp_eta),
    nrow_e_Xq       = NROW(e_mod_stuff$xtemp),
    e_Xq            = e_mod_stuff$xtemp,
    e_times         = c(e_mod_stuff$eventtime, unlist(e_mod_stuff$quadpoint)),
    e_d             = c(e_mod_stuff$d, rep(1, length(unlist(e_mod_stuff$quadpoint)))),
    e_xbar          = as.array(e_mod_stuff$xbar),
    e_weights       = as.array(e_weights),
    e_weights_rep   = as.array(rep(e_weights, times = quadnodes)),
    quadweight      = as.array(e_mod_stuff$quadweight),
    
    # priors
    y_prior_dist = c(unlist(unique(fetch(y_prior_stuff, "prior_dist")))), 
    e_prior_dist = e_prior_stuff$prior_dist,
    a_prior_dist = a_prior_stuff$prior_dist,    
    y_prior_dist_for_intercept = fetch_array(y_prior_intercept_stuff, "prior_dist"),  
    e_prior_dist_for_intercept = e_prior_intercept_stuff$prior_dist,
    y_prior_dist_for_aux = fetch_array(y_prior_aux_stuff, "prior_dist"),  
    e_prior_dist_for_aux = e_prior_aux_stuff$prior_dist,
    
    # hyperparameters for priors
    y_prior_mean = fetch_array(y_prior_stuff, "prior_mean"), 
    e_prior_mean = e_prior_stuff$prior_mean, 
    a_prior_mean = a_prior_stuff$prior_mean, 
    y_prior_mean_for_intercept = fetch_array(y_prior_intercept_stuff, "prior_mean"),
    e_prior_mean_for_intercept = c(e_prior_intercept_stuff$prior_mean),
    y_prior_scale = fetch_array(y_prior_stuff, "prior_scale"), 
    e_prior_scale = e_prior_stuff$prior_scale, 
    a_prior_scale = a_prior_stuff$prior_scale, 
    y_prior_scale_for_intercept = fetch_array(y_prior_intercept_stuff, "prior_scale"), 
    e_prior_scale_for_intercept = c(e_prior_intercept_stuff$prior_scale), 
    y_prior_df = fetch_array(y_prior_stuff, "prior_df"), 
    e_prior_df = e_prior_stuff$prior_df, 
    a_prior_df = a_prior_stuff$prior_df, 
    y_prior_df_for_intercept = fetch_array(y_prior_intercept_stuff, "prior_df"),  
    e_prior_df_for_intercept = c(e_prior_intercept_stuff$prior_df),
    y_global_prior_scale = c(unique(fetch_(y_prior_stuff, "global_prior_scale"))),
    e_global_prior_scale = e_prior_stuff$global_prior_scale,
    a_global_prior_scale = a_prior_stuff$global_prior_scale,
    y_global_prior_df = c(unique(fetch_(y_prior_stuff, "global_prior_df"))), 
    e_global_prior_df = e_prior_stuff$global_prior_df,
    a_global_prior_df = a_prior_stuff$global_prior_df,
    prior_PD = as.integer(prior_PD),
    long_lp = as.integer(long_lp),
    event_lp = as.integer(event_lp)
  )
  
  # data for association structure
  # !!! must be careful with corresponding use of indexing in Stan code
  standata$assoc <- as.integer(a_K > 0L) 

  # indicator for which components are required to build the association terms
  standata$assoc_uses <- sapply(
    c("etavalue", "etaslope", "etalag", "etaauc", "muvalue",  "muslope",  "mulag",  "muauc"), 
    function(x, assoc) {
      nm_check <- switch(x,
                         etavalue = "etavalue|muvalue|etaslope|muslope",
                         etaslope = "etaslope|muslope",
                         etalag   = "etalag|mulag",
                         etaauc   = "etaauc|muauc",
                         muvalue  = "muvalue|muslope",
                         muslope  = "muslope",
                         mulag    = "mulag",
                         muauc    = "muauc")
      sel <- grep(nm_check, rownames(assoc))
      as.integer(any(unlist(assoc[sel,])))
    }, assoc = assoc)  
    
  # indexing
  # 1 = ev; 2 = es; 3 = el; 4 = ea; 5 = mv; 6 = ms; 7 = ml; 8 = ma;
  # 9 = shared_b; 10 = shared_coef;
  # 11 = ev_data; 12 = mv_data; 13 = es_data; 14 = ms_data;
  # 15 = evev; 16 = evmv; 17 = mvev; 18 = mvmv;
  sel <- grep("which|null", rownames(assoc), invert = TRUE)
  standata$has_assoc <- matrix(as.integer(assoc[sel,]), ncol = M) 

  # auc association structure
  standata$auc_quadnodes <- as.integer(auc_quadnodes)
  standata$Npat_times_auc_quadnodes <- as.integer(e_mod_stuff$Npat * auc_quadnodes)  
  standata$nrow_y_Xq_auc <- as.integer(auc_quadnodes * NROW(a_mod_stuff[[1]]$xtemp_eta))
  standata$auc_quadweights <- if (standata$assoc_uses[4]) as.array(auc_quadweights) else double(0)
  
  # shared random effects
  standata$which_b_zindex    <- as.array(unlist(assoc["which_b_zindex",]))
  standata$which_coef_zindex <- as.array(unlist(assoc["which_coef_zindex",]))
  standata$which_coef_xindex <- as.array(unlist(assoc["which_coef_xindex",]))
  standata$size_which_b      <- as.array(sapply(assoc["which_b_zindex",    ], length))
  standata$size_which_coef   <- as.array(sapply(assoc["which_coef_zindex", ], length))
  
  # interactions between association terms
  standata$which_interactions      <- unlist(assoc["which_interactions",])
  standata$size_which_interactions <- c(sapply(assoc["which_interactions",], sapply, length))
 
  # design matrix for interactions between association terms and data
  standata$y_Xq_data <- as.array(t(as.matrix(do.call("cbind", fetch(a_mod_stuff, "xtemp_data")))))
  # number of columns in y_Xq_data corresponding to each interaction type 
  # (ie, etavalue, etaslope, muvalue, muslope) for each submodel
  standata$a_K_data  <- fetch_array(a_mod_stuff, "K_data")  

  # Indexing for combined beta vector, response vector, design matrix, weights, etc
  y_y        <- fetch(y_mod_stuff, "y")
  y_K        <- fetch_(y_mod_stuff, "K")
  y_N        <- fetch_(y_mod_stuff, "N")
  y_real_N   <- fetch_(y_mod_stuff, "real_N")
  y_int_N    <- fetch_(y_mod_stuff, "int_N")
  y_is_real  <- fetch_(y_mod_stuff, "is_real")
  standata$y_real     <- as.array(as.numeric(unlist(y_y[y_is_real])))
  standata$y_int      <- as.array(as.integer(unlist(y_y[!y_is_real])))
  standata$y_K_beg    <- as.array(sapply(1:M, function(m) sum(y_K[0:(m-1)]) + 1))
  standata$y_K_end    <- as.array(sapply(1:M, function(m) sum(y_K[0:m])))
  standata$y_beg      <- as.array(sapply(1:M, function(m) sum(y_N[0:(m-1)]) + 1))
  standata$y_end      <- as.array(sapply(1:M, function(m) sum(y_N[0:m])))
  standata$y_real_beg <- as.array(sapply(1:M, function(m) if (y_is_real[m])  sum(y_real_N[0:(m-1)]) + 1L else 0L))
  standata$y_real_end <- as.array(sapply(1:M, function(m) if (y_is_real[m])  sum(y_real_N[0:m])          else 0L))  
  standata$y_int_beg  <- as.array(sapply(1:M, function(m) if (!y_is_real[m]) sum(y_int_N[0:(m-1)]) + 1L  else 0L))
  standata$y_int_end  <- as.array(sapply(1:M, function(m) if (!y_is_real[m]) sum(y_int_N[0:m])           else 0L))    
    
  # Sum dimensions
  for (i in c("y_K", "y_N", "y_real_N", "y_int_N", "y_has_aux", "a_K_data",
              paste0("y_has_intercept", c("", "_unbound", "_lobound", "_upbound")),
              paste0("size_which_",     c("b", "coef", "interactions")))) {
    standata[[paste0("sum_", i)]] <- as.integer(sum(standata[[i]]))
  }
  
  # product normal
  standata$y_num_normals <-
    if (any(standata$y_prior_dist == 7)) as.integer(standata$y_prior_df) else integer(0)
  standata$e_num_normals <-
    if (standata$e_prior_dist == 7) as.integer(standata$e_prior_df) else integer(0)
  standata$a_num_normals <-
    if (standata$a_prior_dist == 7) as.integer(standata$a_prior_df) else integer(0)
  
  # hyperparameters for priors of auxiliary parameters
  standata$y_prior_mean_for_aux  <- fetch_array(y_prior_aux_stuff, "prior_mean")
  standata$e_prior_mean_for_aux  <- if (basehaz$type == 1L) as.array(0) else as.array(e_prior_aux_stuff$prior_mean)  
  standata$y_prior_scale_for_aux <- fetch_array(y_prior_aux_stuff, "prior_scale")
  standata$e_prior_scale_for_aux <- e_prior_aux_stuff$prior_scale  
  standata$y_prior_df_for_aux    <- fetch_array(y_prior_aux_stuff, "prior_df")
  standata$e_prior_df_for_aux    <- e_prior_aux_stuff$prior_df  
  
  # data for random effects
  group <- lapply(y_mod_stuff, function(x) {
    pad_reTrms(Ztlist = x$Ztlist, 
               cnms   = x$cnms, 
               flist  = x$flist)})
  Z              <- fetch(group, "Z")
  y_cnms         <- fetch(group, "cnms")
  y_flist_padded <- fetch(group, "flist")
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
  standata$t    <- as.integer(t)
  standata$t_i  <- as.integer(t_i)
  standata$pmat <- as.array(pmat)
  standata$p    <- as.array(p)
  standata$l    <- as.array(l)
  standata$qmat <- as.array(qmat)
  standata$q    <- 0L # not used
  standata$q1   <- as.array(q1)
  standata$q2   <- as.array(q2)
  standata$len_theta_L <- sum(choose(p, 2), p)
  Zmerge <- Matrix::bdiag(Z)
  standata$len_b <- as.integer(ncol(Zmerge))
  parts <- rstan::extract_sparse_parts(Zmerge)
  standata$num_non_zero <- as.integer(length(parts$w))
  standata$w <- parts$w
  standata$v <- parts$v
  standata$u <- as.array(parts$u)
  
  # data for calculating eta, eta slope, eta lag, eta auc in GK quadrature 
  for (i in c("eta", "eps", "lag", "auc")) {
    nm_check <- switch(i,
                       eta = "^eta|^mu",
                       eps = "slope",
                       lag = "etalag|mulag",
                       auc = "auc")
    sel <- grep(nm_check, rownames(assoc))
    if (any(unlist(assoc[sel,]))) {
      X_tmp <- as.matrix(Matrix::bdiag(fetch(a_mod_stuff, paste0("xtemp_", i))))
      group_tmp <- lapply(a_mod_stuff, function(x) {
        pad_reTrms(Ztlist = x[[paste0("group_", i)]][["Ztlist"]], 
                   cnms   = x[[paste0("group_", i)]][["cnms"]], 
                   flist  = x[[paste0("group_", i)]][["flist"]])})
      Z_tmp <- Matrix::bdiag(fetch(group_tmp, "Z"))      
    } else {
      X_tmp <- matrix(0,0,standata$sum_y_K)
      Z_tmp <- matrix(0,0,0) 
    }
    parts_Z_tmp <- rstan::extract_sparse_parts(Z_tmp)
    standata[[paste0("y_Xq_", i)]] <- as.array(X_tmp)
    standata[[paste0("num_non_zero_Zq_", i)]] <- as.integer(length(parts_Z_tmp$w))
    standata[[paste0("w_Zq_", i)]] <- parts_Z_tmp$w
    standata[[paste0("v_Zq_", i)]] <- parts_Z_tmp$v
    standata[[paste0("u_Zq_", i)]] <- as.array(parts_Z_tmp$u)    
  }
  
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
    if (y_mod_stuff[[x]]$is_bernoulli) return_fam <- 4L
    return_fam}))
  standata$any_fam_3 <- as.integer(any(standata$family == 3L))
  
  # B-splines baseline hazard
  if (basehaz$type_name == "weibull") {
    standata$e_basehaz_X <- matrix(log(standata$e_times), length(standata$e_times), 1) 
  } else if (basehaz$type_name == "bs") {
    standata$e_basehaz_X <- as.array(predict(basehaz$bs_basis, standata$e_times)) 
  } else if (basehaz$type_name == "piecewise") {
    e_times_quantiles <- cut(standata$e_times, basehaz$knots, include.lowest = TRUE, labels = FALSE)
    tmp <- matrix(NA, length(e_times_quantiles), basehaz$df)
    for (i in 1:basehaz$df) tmp[, i] <- ifelse(e_times_quantiles == i, 1, 0)
    standata$e_basehaz_X <- as.array(tmp)
  } else {
    standata$e_basehaz_X <- matrix(0,0,0)  
  }
  
  #================
  # Initial values
  #================
  
  if (init == "model_based")
    init <- generate_init_function(y_mod_stuff, e_mod_stuff, standata)

  #===========
  # Fit model
  #===========
    
  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$jm
  pars <- c(if (standata$sum_y_has_intercept_unbound) "y_gamma_unbound",
            if (standata$sum_y_has_intercept_lobound) "y_gamma_lobound",
            if (standata$sum_y_has_intercept_upbound) "y_gamma_upbound", 
            if (standata$sum_y_K) "y_beta",
            if (standata$e_has_intercept) "e_gamma",
            if (standata$e_K) "e_beta",
            if (standata$a_K) "a_beta",
            if (standata$Npat) "b_by_model",
            if (standata$sum_y_has_aux) "y_aux",
            "e_aux",
            if (standata$len_theta_L) "theta_L")
            
  cat(paste0(if (M == 1L) "Uni" else "Multi", "variate joint model specified\n"))
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
  
  # Names for pars
  y_nms <- unlist(lapply(1:M, function(m) 
    if (ncol(y_mod_stuff[[m]]$xtemp)) paste0("Long", m, "|", colnames(y_mod_stuff[[m]]$xtemp))))
  y_int_nms <- unlist(lapply(1:M, function(m) 
    if (y_mod_stuff[[m]]$has_intercept) paste0("Long", m, "|(Intercept)")))
  y_aux_nms <- character()  
  for (m in 1:M) {
    if (is.gaussian(y_mod_stuff[[m]]$famname))   y_aux_nms <- c(y_aux_nms, paste0("Long", m,"|sigma"))
    else if (is.gamma(y_mod_stuff[[m]]$famname)) y_aux_nms <- c(y_aux_nms, paste0("Long", m,"|shape"))
    else if (is.ig(y_mod_stuff[[m]]$famname))    y_aux_nms <- c(y_aux_nms, paste0("Long", m,"|lambda"))
    else if (is.nb(y_mod_stuff[[m]]$famname))    y_aux_nms <- c(y_aux_nms, paste0("Long", m,"|reciprocal_dispersion"))
  }  
  e_nms <- if (ncol(e_mod_stuff$xtemp)) paste0("Event|", colnames(e_mod_stuff$xtemp))    
  e_int_nms <- if (e_has_intercept) "Event|(Intercept)"
  e_aux_nms <- if (basehaz$type == 1L) "Event|weibull-shape" else paste0("Event|basehaz-coef", seq(basehaz$df))               
  a_nms <- character()  
  for (m in 1:M) {
    if (assoc["etavalue",         ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etavalue"))
    if (assoc["etavalue_data",    ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etavalue:", colnames(a_mod_stuff[[m]][["xq_data"]][["etavalue_data"]])))
    if (assoc["etavalue_etavalue",][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etavalue:Long", assoc["which_interactions",][[m]][["etavalue_etavalue"]], "|etavalue"))
    if (assoc["etavalue_muvalue", ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etavalue:Long", assoc["which_interactions",][[m]][["etavalue_muvalue"]],  "|muvalue"))
    if (assoc["etaslope",         ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etaslope"))
    if (assoc["etaslope_data",    ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etaslope:", colnames(a_mod_stuff[[m]][["xq_data"]][["etaslope_data"]])))    
    if (assoc["etalag",           ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etalag"))
    if (assoc["etaauc",           ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etaauc"))
    if (assoc["muvalue",          ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue"))
    if (assoc["muvalue_data",     ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue:", colnames(a_mod_stuff[[m]][["xq_data"]][["muvalue_data"]])))    
    if (assoc["muvalue_etavalue", ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_etavalue"]], "|etavalue"))
    if (assoc["muvalue_muvalue",  ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_muvalue"]],  "|muvalue"))
    if (assoc["muslope",          ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muslope"))
    if (assoc["muslope_data",     ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muslope:", colnames(a_mod_stuff[[m]][["xq_data"]][["muslope_data"]])))    
    if (assoc["mulag",            ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|mulag"))
    if (assoc["muauc",            ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muauc"))
  }
  if (sum(standata$size_which_b)) {
    temp_g_nms <- lapply(1:M, FUN = function(m) {
      all_nms <- paste0(paste0("Long", m, "|b["), y_mod_stuff[[m]]$cnms[[id_var]], "]")
      all_nms[assoc["which_b_zindex",][[m]]]})
    a_nms <- c(a_nms, paste0("Assoc|", unlist(temp_g_nms)))
  }
  if (sum(standata$size_which_coef)) {
    temp_g_nms <- lapply(1:M, FUN = function(m) {
      all_nms <- paste0(paste0("Long", m, "|coef["), y_mod_stuff[[m]]$cnms[[id_var]], "]")
      all_nms[assoc["which_coef_zindex",][[m]]]})
    a_nms <- c(a_nms, paste0("Assoc|", unlist(temp_g_nms)))
  }
  
  # Sigma names
  if (standata$len_theta_L) {
    thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                      permuted = FALSE)
    nc <- sapply(cnms, FUN = length)
    nms <- names(cnms)
    Sigma <- apply(thetas, 1:2, FUN = function(theta) {
      Sigma <- mkVarCorr(sc = 1, cnms, nc, theta, nms)
      unlist(sapply(Sigma, simplify = FALSE, 
                    FUN = function(x) x[lower.tri(x, TRUE)]))
    })
    l <- length(dim(Sigma))
    end <- tail(dim(Sigma), 1L)
    shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
    if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
      stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain] 
    }
    else for (chain in 1:end) {
      stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
    }
    Sigma_nms <- lapply(cnms, FUN = function(grp) {
      nm <- outer(grp, grp, FUN = paste, sep = ",")
      nm[lower.tri(nm, diag = TRUE)]
    })
    for (j in seq_along(Sigma_nms)) {
      Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
    }
    Sigma_nms <- unlist(Sigma_nms)
  }
  
  
  new_names <- c(y_int_nms,
                 y_nms,
                 e_int_nms,
                 e_nms,
                 a_nms,                   
                 if (length(group)) c(paste0("b[", b_nms, "]")),
                 y_aux_nms,
                 e_aux_nms,
                 if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                 "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  
  n_grps <- standata$l - 1
  names(n_grps) <- cnms_nms  # n_grps is num. of levels within each grouping factor
  names(p) <- cnms_nms       # p is num. of variables within each grouping factor
  
  # Undo ordering of matrices if bernoulli
  y_mod_stuff <- lapply(y_mod_stuff, unorder_bernoulli)
  
  fit <- nlist(stanfit, family, formula = c(formulaLong, formulaEvent), 
               id_var, time_var, offset = NULL, weights, quadnodes, basehaz,
               M, cnms, n_yobs = unlist(list_nms(fetch(y_mod_stuff, "N"), M)), 
               n_subjects = e_mod_stuff$Npat, n_grps, assoc,
               y_mod_stuff, e_mod_stuff, fr = c(fetch(a_mod_stuff, "fr"), list(e_mod_stuff$model_frame)),
               y = list_nms(fetch(y_mod_stuff, "y"), M),
               d = e_mod_stuff$d, eventtime = e_mod_stuff$eventtime,
               epsilon = if (standata$assoc_uses[2]) eps else NULL,
               dataLong, dataEvent, call, na.action, algorithm, 
               glmod = fetch(y_mod_stuff, "mod"), coxmod = e_mod_stuff$mod,
               standata = NULL, terms = NULL, prior.info = NULL)  # last line elements not currently used
  out <- stanjm(fit)
  
  return(out)
}


#--------------- Functions related to longitudinal submodel

# Fit separate longitudinal submodel
#
# @param mc The (slightly modified) user specified call for the longitudinal submodel
# @param family The family object for the longitudinal submodel
# @param supported_families A string vector of supported family names
# @param supported_links A string vector of supported link names
# @param sparse A logical specifying whether to use a sparse design matrix
# @return A named list
handle_glmod <- function(mc, family, supported_families, supported_links, sparse) {
  
  if ((family$family == "gaussian") && (family$link == "identity")) {
    mc[[1]]    <- quote(lme4::lmer)
    mc$family  <- NULL
    mc$control <- get_control_args()
  } else if (family$family == "neg_binomial_2") {
    mc[[1]]    <- quote(lme4::glmer.nb)
    mc$family  <- NULL
    mc$control <- get_control_args(glmer = TRUE)               
  } else {
    mc[[1]]    <- quote(lme4::glmer)
    mc$control <- get_control_args(glmer = TRUE)               
  }
  mod <- eval(mc, parent.frame())      	
  
  # Response vector
  y <- as.vector(lme4::getME(mod, "y"))
  y <- validate_glm_outcome_support(y, family)
  if (is.binomial(family$family) && NCOL(y) == 2L) {
    y      <- as.integer(y[, 1L])
    trials <- as.integer(y[, 1L] + y[, 2L])
  } else {
    trials <- rep(0L, length(y))
  }
  
  # Design matrix
  x <- as.matrix(lme4::getME(mod, "X"))
  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  if (has_intercept) {
    check_int <- check_intercept(family$family, family$link)
    has_intercept_unbound <- check_int$unbound
    has_intercept_lobound <- check_int$lobound
    has_intercept_upbound <- check_int$upbound
  } else {
    has_intercept_unbound <- 0L
    has_intercept_lobound <- 0L
    has_intercept_upbound <- 0L
  } 
  offset <- lme4::getME(mod, "offset")
  
  # Random effect terms
  group_call <- mc
  group_call[[1]] <- quote(lme4::glFormula)
  group <- eval(group_call, parent.frame())$reTrms      	
  
  Ztlist <- group$Ztlist
  cnms   <- group$cnms
  flist  <- group$flist
  
  # Reorder y, X, Z if bernoulli (zeros first)
  if (is.binomial(family$family) && all(y %in% 0:1)) {      
    ord    <- order(y)
    y      <- y     [ord]
    trials <- trials[ord]
    xtemp  <- xtemp [ord, , drop = FALSE]  
    Ztlist <- lapply(Ztlist, function(x) x[, ord, drop = FALSE]) 
    flist  <- lapply(flist,  function(x) x[ord]) 
    N01    <- sapply(0:1,    function(x) sum(y == x))
  } else {
    ord    <- NULL
    N01    <- rep(0L, 2)  # dud entry if not bernoulli
  }
  
  # Dimensions
  N       <- NROW(xtemp)
  is_real <- check_response_real(family$family)
  real_N  <- if (is_real) N else 0L
  int_N   <- if (!is_real) N else 0L
  K       <- ncol(xtemp)	
  has_aux <- check_for_aux(family$family)
  
  # Update formula if using splines or other data dependent predictors
  vars <- get_formvars(mod)

  # Family indicators
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  famname <- supported_families[fam]
  is_bernoulli <- is.binomial(famname) && all(y %in% 0:1)
  is_nb <- is.nb(famname)
  is_gaussian <- is.gaussian(famname)
  is_gamma <- is.gamma(famname)
  is_ig <- is.ig(famname)
  is_continuous <- is_gaussian || is_gamma || is_ig
  is_lmer <- is.lmer(family)
  
  # Require intercept for certain family and link combinations
  if (!has_intercept) {
    linkname <- supported_links[link]
    needs_intercept <- !is_gaussian && linkname == "identity" ||
      is_gamma && linkname == "inverse" ||
      is.binomial(famname) && linkname == "log"
    if (needs_intercept)
      stop("To use the specified combination of family and link (", famname, 
           ", ", linkname, ") the model must have an intercept.")
  }    
    
  # Return list
  nlist(mod, vars, is_real, y, x, xtemp, xbar, trials, N01, ord, famname,
    offset, Ztlist, cnms, flist, has_intercept, has_intercept_unbound,
    has_intercept_lobound, has_intercept_upbound, has_aux, N, real_N, int_N, K,
    is_bernoulli, is_nb, is_gaussian, is_gamma, is_ig, is_continuous, is_lmer)
}

# Function to return a single cnms object for all longitudinal submodels
#
# @param x A list, with each element being a cnms object returned by (g)lmer
get_common_cnms <- function(x) {
  nms <- lapply(x, names)
  unique_nms <- unique(unlist(nms))
  cnms <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    unlist(lapply(1:length(x), function(m) 
      if (nm %in% nms[[m]]) paste0("Long", m, "|", x[[m]][[nm]])))
  })
  names(cnms) <- unique_nms
  cnms
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
# @param mod A (g)lmer model object from which to extract the formula and terms
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
# @param id_var The name of the ID variable
# @param y_flist A list containing the flist objects returned for each 
#   separate longitudinal submodel
# @return A vector of factor levels corresponding to the IDs appearing
#   in the longitudinal submodels
check_id_list <- function(id_var, y_flist) {
  id_list <- unique(lapply(y_flist, function(x) levels(x[[id_var]])))
  if (length(id_list) > 1L)
    stop("The subject IDs are not the same in all longitudinal submodels.",
         call. = FALSE)
  unlist(id_list)  
}

# Function to return control arguments for lmer or glmer call
#
# @param glmer A logical indicating whether to use glmerControl
# @param norank
get_control_args <- function(glmer = FALSE, norank = FALSE) {
  do.call(ifelse(glmer, "glmerControl", "lmerControl"),
    list(
      check.nlev.gtreq.5 = "ignore",
      check.nlev.gtr.1 = "stop",
      check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nlev = "ignore",
      check.nobs.vs.nRE = "ignore",
      check.rankX = ifelse(norank, "ignore", "message+drop.cols")))
}

# Undo the ordering of the response vector, design matrices, weights, etc
# if they were a bernoulli response sorted for model estimation
#
# @param mod_stuff A named list returned by a call to handle_glmod
unorder_bernoulli <- function(mod_stuff) {
  if (mod_stuff$is_bernoulli) {
    if (is.null(mod_stuff$ord))
      stop("Bernoulli sorting vector not found.")
    mod_stuff$y       <- mod_stuff$y[order(mod_stuff$ord)]
    mod_stuff$trials  <- mod_stuff$trials[order(mod_stuff$ord)]
    mod_stuff$weights <- mod_stuff$weights[order(mod_stuff$ord)]
    mod_stuff$xtemp   <- mod_stuff$xtemp[order(mod_stuff$ord), , drop = FALSE]
    mod_stuff$Ztlist  <- lapply(mod_stuff$Ztlist, function(x) x[, order(mod_stuff$ord), drop = FALSE])
    mod_stuff$flist   <- lapply(mod_stuff$flist,  function(x) x[order(mod_stuff$ord)])
  }
  mod_stuff
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

# Function to check if the submodel should include a auxiliary term
#
# @param family A GLM family
# @return A logical specify whether the submodel includes a auxiliary term
check_for_aux <- function(family) {
  if (family == "binomial" || family == "poisson") FALSE else TRUE
}


#--------------- Functions related to event submodel

# Fit separate event submodel
#
# @param mc The (slightly modified) user specified call for the event submodel
# @param quadnodes An integer, the user specified number of GK quadrature nodes
# @param id_var The name of the ID variable
# @param unique_id_list A character vector with the unique IDs (factor levels)
#   that appeared in the longitudinal submodels
# @param sparse A logical indicating whether to use a sparse design matrix
handle_coxmod <- function(mc, quadnodes, id_var, unique_id_list, sparse) {
  
  mc[[1]] <- quote(survival::coxph) 
  mc$x    <- TRUE
  mod <- eval(mc, parent.frame())
  mf1 <- mf2 <- expand.model.frame(mod, id_var, na.expand = TRUE)
  mf2 <- cbind(unclass(mf2[,1]), mf2[, -1, drop = FALSE])
  y   <- mod$y
  
  # Entry and exit times
  entrytime <- rep(0, length(unique_id_list)) # entry times currently assumed to be zero for all individuals, i.e. no delayed entry
  if (attr(y, "type") == "counting") {
    tvc         <- TRUE
    mf_event    <- do.call(rbind, lapply(split(mf2, mf2[, id_var]), function(d) d[which.max(d[,"stop"]), ]))
    eventtime   <- mf_event[["stop"]]
  } else if (attr(y, "type") == "right") {
    tvc         <- FALSE 
    mf_event    <- mf2
    eventtime   <- mf_event[["time"]]
  } else stop("Only 'right' or 'counting' type Surv objects are allowed 
               on the LHS of the event submodel formula")
  
  # Event indicator and ID list
  d     <- mf_event[["status"]]  
  flist <- mf_event[[id_var]]
  names(eventtime) <- names(d) <- flist
  
  # Error checks for the ID variable
  if (!identical(unique_id_list, levels(factor(flist))))
    stop("The patient IDs (levels of the grouping factor) included ",
         "in the longitudinal and event submodels do not match")
  if (is.unsorted(factor(flist)))
    stop("'dataEvent' needs to be sorted by the subject ",
         "ID/grouping variable", call. = FALSE)
  if (!identical(length(unique_id_list), length(eventtime)))
    stop("The number of patients differs between the longitudinal and ",
         "event submodels. Perhaps you intended to use 'start/stop' notation ",
         "for the Surv() object.")
  
  # Unstandardised quadrature points
  quadpoint <- lapply(get_quadpoints(quadnodes)$points, unstandardise_quadpoints, entrytime, eventtime)
  quadtimes <- c(list(eventtime), quadpoint)
  names(quadpoint) <- paste0("quadpoint", seq(quadnodes))
  names(quadtimes) <- c("eventtime", names(quadpoint))
  
  # Obtain design matrix at event times and unstandardised quadrature points
  
  if (tvc) {  # time varying covariates in event model
    
    # Model frame at event times
    mf2           <- data.table(mf2, key = c(id_var, "start"))
    mf2[["start"]] <- as.numeric(mf2[["start"]])
    mf2_eventtime <- mf2[, .SD[.N], by = get(id_var)]
    mf2_eventtime <- mf2_eventtime[, get := NULL]   
    
    # Model frame corresponding to observation times which are 
    #   as close as possible to the unstandardised quadrature points                      
    mf2_q  <- do.call(rbind, lapply(quadpoint, FUN = function(x)
      mf2[data.table::SJ(flist, x), roll = TRUE, rollends = c(TRUE, TRUE)]))
    
    # Model frame evaluated at both event times and quadrature points
    mf2_q <- rbind(mf2_eventtime, mf2_q)
    
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are time varying covariates in the event submodel and
    #   therefore the design matrix differs depending on the quadrature point 
    fm_RHS <- delete.response(terms(mod))
    x_quadtime   <- model.matrix(fm_RHS, data = mf2_q)
    
  } else {  # no time varying covariates in event model
    
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are no time varying covariates in the event submodel and
    #   therefore the design matrix is identical at event time and at all
    #   quadrature points
    x_quadtime   <- do.call(rbind, lapply(1:(quadnodes + 1), FUN = function(x) mod$x))
    
  }
  
  # Centering of design matrix for event model
  x <- as.matrix(x_quadtime) 
  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  
  # Some additional bits -- NB weights here are quadrature weights, not prior weights
  K <- NCOL(xtemp)
  Npat <- length(eventtime)
  quadweight <- unlist(lapply(get_quadpoints(quadnodes)$weights, unstandardise_quadweights, entrytime, eventtime))
  
  # Return list
  nlist(mod, entrytime, eventtime, d, x, xtemp, xbar, flist, 
        quadpoint, quadweight, quadtimes, tvc, K, Npat, model_frame = mf1)
}

# Deal with the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements df, knots 
# @param ok_basehaz A list of admissible baseline hazards
# @param eventtime A numeric vector with eventtimes for each individual
# @param d A numeric vector with event indicators for each individual
handle_basehaz <- function(basehaz, basehaz_ops, 
                           ok_basehaz = nlist("weibull", "bs", "piecewise"),
                           ok_basehaz_ops = nlist("df", "knots"),
                           eventtime, d) {

  if (!basehaz %in% unlist(ok_basehaz))
    stop("The baseline hazard should be one of ", paste(names(ok_basehaz), collapse = ", "))
  if (!all(names(basehaz_ops) %in% unlist(ok_basehaz_ops)))
    stop("The baseline hazard options list can only include ", paste(names(ok_basehaz_ops), collapse = ", "))
  
  type      <- switch(basehaz,
                      weibull = 1L,
                      bs = 2L,
                      piecewise = 3L)
  type_name <- basehaz
  user_df   <- basehaz_ops$df
  df        <- basehaz_ops$df
  knots     <- basehaz_ops$knots
  bs_basis  <- NULL
  
  if (type_name == "weibull") {
    # handle df and knots
    if (!is.null(df))
      warning("'df' will be ignored since baseline hazard was set to weibull.", 
              immediate. = TRUE, call. = FALSE)
    if (!is.null(knots))
      warning("'knots' will be ignored since baseline hazard was set to weibull.", 
              immediate. = TRUE, call. = FALSE) 
    user_df <- NULL
    df      <- 1L
    knots   <- NULL
  } else if (type_name %in% c("bs", "piecewise")) {
    # handle df and knots
    if (!any(is.null(df), is.null(knots))) { 
      # both specified
      stop("Cannot specify both 'df' and 'knots' for the baseline hazard.", call. = FALSE)
    } else if (all(is.null(df), is.null(knots))) { 
      # both null -- use default df
      user_df <- df <- 6L
      knots <- NULL
    } else if (!is.null(df)) { 
      # only df specified
      if (type == 2L) {
        if (df < 3) stop("'df' must be at least 3 for B-splines baseline hazard.")
        user_df <- df <- df + 1
      }
    } else if (!is.null(knots)) {          
      # only knots specified
      if (!is.numeric(knots)) stop("'knots' vector must be numeric", call. = FALSE)
      if (any(knots < 0)) stop("'knots' must be non-negative", call. = FALSE)      
      if (type == 2L) df <- length(knots) + 4
      else if (type == 3L) df <- length(knots) + 1
    } else {
      stop("Bug found: unable to reconcile 'df' and 'knots' arguments.", call. = FALSE) 
    }
  }  
 
  # Evaluate spline basis (knots, df, etc) based on distribution of observed event times
  # or evaluate cut points for piecewise constant baseline hazard
  if (type == 2L) {
    bs_basis <- splines::bs(eventtime[(d > 0)], df = user_df, knots = knots, 
                            Boundary.knots = c(0, max(eventtime)), intercept = TRUE)
  } else if (type == 3L) {
    if (is.null(knots)) {
      knots <- quantile(eventtime[(d > 0)], probs = seq(0, 1, 1 / df))
      knots[[1]] <- 0
      knots[[length(knots)]] <- max(eventtime)
    } else {
      if (any(knots > max(eventtime)))
        stop("'knots' for the baseline hazard cannot be greater than the ",
             "largest event time.", call. = FALSE)
      knots <- c(0, knots, max(eventtime))
    }
  }  
   
  nlist(type, type_name, user_df, df, knots, bs_basis)   
}

# Function to return standardised GK quadrature points and weights
#
# @param nodes The required number of quadrature nodes
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


#--------------- Functions related to association structure

# Return a named list with information about the specified association structure 
# 
# @param user_x A character vector or NULL, being the user input to the
#   assoc argument (for one submodel) in the stan_jm call
# @param y_mod_stuff A list returned by a call to handle_glmod
# @param id_var The name of the ID variable 
# @param M Integer specifying the total number of longitudinal submodels
# @return A list with information about the desired association structure
validate_assoc <- function(user_x, y_mod_stuff, ok_assoc, ok_assoc_data,
                           ok_assoc_interactions, id_var, M) {

  ok_inputs <- c(ok_assoc, paste0(ok_assoc_data, "_data"),
                 unlist(lapply(ok_assoc_interactions, paste0, "_", ok_assoc_interactions))) 

  # Check user input to assoc argument
  trimmed_x <- trim_assoc(user_x, ok_assoc_data, ok_assoc_interactions)
  if (is.null(user_x) || all(trimmed_x %in% ok_inputs)) {
    assoc <- sapply(ok_inputs, `%in%`, trimmed_x, simplify = FALSE)
    if (is.null(user_x)) {
      assoc$null <- TRUE
    } else if (is.vector(user_x) && is.character(user_x)) {
      if ((assoc$null) && (length(user_x) > 1L))
        stop("In assoc, 'null' cannot be specified in conjuction ",
             "with another association type", call. = FALSE)
      STOP_combination_not_allowed(assoc, "etavalue", "muvalue")
      STOP_combination_not_allowed(assoc, "etaslope", "muslope")
      STOP_combination_not_allowed(assoc, "etalag",   "mulag")
    } else {
      stop("'assoc' argument should be a character vector or, for a multivariate ",
           "joint model, possibly a list of character vectors.", call. = FALSE)    
    }    
  } else {
    stop("An unsupported association type has been specified. The ",
         "'assoc' argument can only include the following association ", 
         "types: ", paste(ok_assoc, collapse = ", "), ", as well as ",
         "possible interactions either between association terms or ",
         "with observed data.", call. = FALSE)  
  }
  
  # Parse suffix specifying desired time lag
  ok_inputs_lag <- c("etalag", "mulag")
  assoc$which_lag <- parse_assoc_lag(ok_inputs_lag, user_x)
       
  # Parse suffix specifying indices for shared random effects
  cnms <- y_mod_stuff$cnms[[id_var]] # names of random effect terms
  assoc$which_b_zindex    <- parse_assoc_sharedRE("shared_b",    user_x, max_index = length(cnms), cnms)
  assoc$which_coef_zindex <- parse_assoc_sharedRE("shared_coef", user_x, max_index = length(cnms), cnms)
  
  if (length(intersect(assoc$which_b_zindex, assoc$which_coef_zindex)))
    stop("The same random effects indices should not be specified in both ",
         "'shared_b' and 'shared_coef'. Specifying indices in 'shared_coef' ",
         "will include both the fixed and random components.", call. = FALSE)
  
  if (length(assoc$which_coef_zindex)) {
    if (length(y_mod_stuff$cnms) > 1L)
      stop("'shared_coef' association structure cannot be used when there is ",
           "clustering at levels other than the individual-level.", call. = FALSE)
    b_nms <- names(assoc$which_coef_zindex)
    assoc$which_coef_xindex <- sapply(b_nms, function(y, beta_nms) {
      beta_match <- grep(y, beta_nms, fixed = TRUE)
      if (!length(beta_match)) {
        stop("In association structure 'shared_coef', no matching fixed effect ",
             "component could be found for the following random effect: ", y, 
             ". Perhaps consider using 'shared_b' association structure instead.")
      } else if (length(beta_match) > 1L) {
        stop("Bug found: In association structure 'shared_coef', multiple ",
             "fixed effect components have been found to match the following ",
             "random effect: ", y)
      }  
      beta_match
    }, beta_nms = colnames(y_mod_stuff$x))
  } else assoc$which_coef_xindex <- numeric(0)
  
  if (!identical(length(assoc$which_coef_zindex), length(assoc$which_coef_xindex)))
    stop("Bug found: the lengths of the fixed and random components of the ",
         "'shared_coef' association structure are not the same.")

  # Parse suffix specifying formula for interactions with data
  ok_inputs_data <- paste0(ok_assoc_data, "_data")
  assoc$which_formulas <- sapply(ok_inputs_data, parse_assoc_data, user_x, simplify = FALSE) 
  
  # Parse suffix specifying indices for interactions between association terms
  ok_inputs_interactions <- unlist(lapply(ok_assoc_interactions, paste0, "_", ok_assoc_interactions))
  assoc$which_interactions <- sapply(ok_inputs_interactions, parse_assoc_interactions, 
                                     user_x, max_index = M, simplify = FALSE)

  assoc
}

# Remove suffixes from the user inputted assoc argument
#
# @param x A character vector, being the user input to the 
#   assoc argument in the stan_jm call
# @param ok_assoc_data A character vector specifying which types
#   of association terms are allowed to be interacted with data
# @param ok_assoc_interactions A character vector specifying which types
#   of association terms are allowed to be interacted with other 
#   association terms
trim_assoc <- function(x, ok_assoc_data, ok_assoc_interactions) {
  x <- gsub("^etalag\\(.*",      "etalag",      x) 
  x <- gsub("^mulag\\(.*",       "mulag",       x) 
  x <- gsub("^shared_b\\(.*",    "shared_b",    x) 
  x <- gsub("^shared_coef\\(.*", "shared_coef", x) 
  for (i in ok_assoc_data)
    x <- gsub(paste0("^", i, "_data\\(.*"),    paste0(i, "_data"), x)
  for (i in ok_assoc_interactions) for (j in ok_assoc_interactions)
    x <- gsub(paste0("^", i, "_", j, "\\(.*"), paste0(i, "_", j),  x) 
  x     
}

# Parse the suffix specified for lagged association terms
# 
# @param ok_inputs A character vector specifying which association
#   structures are allowed to include a suffix specifying a lag time
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @return A numeric vector of length one specifying the time lag to use
parse_assoc_lag <- function(ok_inputs, user_x) {
  val <- grep(paste0("^", ok_inputs, ".*", collapse = "|"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, paste0(ok_inputs, collapse = "|")))[-1]
    if (length(val2) == 1L) {
      lag <- tryCatch(eval(parse(text = paste0("c", val2))), error = function(e) 
        stop("Lagged association structure was specified incorrectly. It should ",
             "include a suffix with the desired lag inside parentheses. See ",
             "Examples in the help file.", call. = FALSE))
      if (length(lag) > 1L) 
        stop("Currently only one lag time is allowed for the lagged association ",
             "structure.", call. = FALSE)
      return(lag)
    } else if (length(val2) > 1L) {
      stop("Only one lagged association structure can be specified for each ",
           "longitudinal submodel.", call. = FALSE)
    } else {
      stop("Lagged association structure was specified incorrectly. It should ",
           "include a suffix with the desired lag inside parentheses. See the ",
           "help file for details.", call. = FALSE)     
    }
  } else numeric(0) 
}

# Parse the formula for specifying a data interaction with an association term
#
# @param x A character string corresponding to one of the allowed
#   association structures for interactions with data, for example, 
#   "etavalue_data" or "etaslope_data"
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @return The parsed formula (which can be used for constructing a 
#   design matrix for interacting data with association type x) or NULL
parse_assoc_data <- function(x, user_x) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    fm <- tryCatch(eval(parse(text = val2)), error = function(e) 
      stop(paste0("Incorrect specification of the formula in the '", x,
                  "' association structure. See Examples in the help file."), call. = FALSE))
    if (!is(fm, "formula"))
      stop(paste0("Suffix to '", x, "' association structure should include ",
                  "a formula within parentheses."), call. = FALSE)
    if (identical(length(fm), 3L))
      stop(paste0("Formula specified for '", x, "' association structure should not ",
                  "include a response."), call. = FALSE)
    if (length(lme4::findbars(fm)))
      stop(paste0("Formula specified for '", x, "' association structure should only ",
                  "include fixed effects."), call. = FALSE)
    if (fm[[2L]] == 1)
      stop(paste0("Formula specified for '", x, "' association structure cannot ",
                  "be an intercept only."), call. = FALSE)
    return(fm)
  } else numeric(0)
}

# Parse the indices specified for shared random effects
#
# @param x A character string corresponding to one of the allowed
#   association structures for shared random effects
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @param max_index An integer specifying the total number of random effects
#   in the longitudinal submodel, and therefore the maximum allowed index for
#   the shared random effects
# @param cnms The names of the random effects corresponding to the 
#   individual-level (id_var) of clustering
# @return A numeric vector specifying indices for the shared random effects
parse_assoc_sharedRE <- function(x, user_x, max_index, cnms) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    if (length(val2)) {
      index <- tryCatch(eval(parse(text = paste0("c", val2))), error = function(e) 
        stop("Incorrect specification of the '", x, "' association structure. ",
             "See Examples in help file.", call. = FALSE))
      if (any(index > max_index))
        stop(paste0("The indices specified for the '", x, "' association structure are ",
                    "greater than the number of subject-specific random effects."), call. = FALSE)
    } else index <- seq_len(max_index)
    names(index) <- cnms[index]
    return(index)   
  } else numeric(0)
}

# Parse the indices specified for interactions between association terms
#
# @param x A character string corresponding to one of the allowed
#   association structures
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @param max_index An integer specifying the maximum allowed index
# @return A numeric vector specifying indices
parse_assoc_interactions <- function(x, user_x, max_index) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    if (length(val2)) {
      index <- tryCatch(eval(parse(text = paste0("c", val2))), error = function(e) 
        stop("Incorrect specification of the '", x, "' association structure. It should ",
             "include a suffix with parentheses specifying the indices of the association ",
             "terms you want to include in the interaction. See Examples in the help file.", call. = FALSE))
      if (any(index > max_index))
        stop("The indices specified for the '", x, "' association structure ",
             "cannot be greater than the number of longitudinal submodels.", call. = FALSE)     
      return(index)
    } else
      stop("Incorrect specification of the '", x, "' association structure. It should ",
             "include a suffix with parentheses specifying the indices of the association ",
             "terms you want to include in the interaction. See Examples in the help file.", call. = FALSE)
  } else numeric(0)      
}

# Make sure that interactions between association terms (for example
# etavalue_etaslope or mu_value_muvalue etc) are always ordered so that
# the first listed association term is for the submodel with the smallest
# index. For example, etavalue1_etavalue2 NOT etavalue2_etavalue1. This
# is to ensure there is no replication such as including both 
# etavalue1_etavalue2 AND etavalue2_etavalue1 when passing to Stan.
#
# @param assoc A named list of named lists, returned by a call to validate_assoc
# @param ok_assoc_interactions A character vector, specifying which association
#   structures are allowed to be used in interactions
check_order_of_assoc_interactions <- function(assoc, ok_assoc_interactions) {
  M <- ncol(assoc)
  for (i in ok_assoc_interactions) {
    for (j in ok_assoc_interactions) {
      header <- paste0(i, "_", j)
      header_reversed <- paste0(j, "_", i)
      for (m in 1:M) {
        if (assoc[header,][[m]]) {
          indices <- assoc["which_interactions",][[m]][[header]]
          sel <- which(indices < m)
          if (length(sel)) {
            # Remove indices for submodels before the current submodel m
            new_indices <- indices[-sel]
            assoc["which_interactions", ][[m]][[header]] <- new_indices
            assoc[header,][[m]] <- (length(new_indices) > 0L)
            # Replace those indices by reversing the order of association terms
            for (k in indices[sel]) {
              assoc["which_interactions",][[k]][[header_reversed]] <- 
                unique(c(assoc["which_interactions",][[k]][[header_reversed]], m))
              assoc[header_reversed,][[k]] <- 
                (length(assoc["which_interactions",][[k]][[header_reversed]]) > 0L)
            }
          }
        }
      }       
    }
  }
  assoc
}

# Return design matrices for evaluating longitudinal submodel quantities 
# at specified quadrature points/times
#
# @param m Integer specifying the longitudinal submodel that the association
#   structure is related to
# @param mc The matched call for the longitudinal submodel
# @param y_mod_stuff A named list returned by a call to handle_glmod (the
#   fit for a single longitudinal submodel)
# @param assoc A named list returned by a call to validate_assoc (details
#   on the desired association structure for a single longitudinal submodel)
# @param has_assoc A list with with indicators for the desired association 
#   structures across all longitudinal submodels 
# @param id_var The name on the ID variable
# @param time_var The name of the time variable
# @param eps The time shift used for the numerical derivative calculation 
#   based on a one-sided different
handle_a_mod <- function(m, mc, y_mod_stuff, e_mod_stuff, assoc, 
                         id_var, time_var, eps) {
  
  e_flist <- e_mod_stuff$flist
  quadtimes <- e_mod_stuff$quadtimes
  
  # Update longitudinal call to include time variable in formula
  mc_stored  <- mc
  mc[[1]]    <- quote(lme4::glFormula)
  mc$formula <- do.call(update.formula, list(mc$formula, paste0("~ . +", time_var)))
  mc$control <- get_control_args(glmer = !y_mod_stuff$is_lmer, norank = TRUE)
  
  # Obtain a model frame with time variable definitely included
  fr <- eval(mc, parent.frame())$fr
  mf <- data.table::data.table(fr, key = c(id_var, time_var)) # create data.table
  mf[[time_var]] <- as.numeric(mf[[time_var]]) # ensure no rounding on merge
  
  # Revert to original formula, and update to reflect predvars (for 
  # example if an individual used data depend predictor variables such
  # as splines with knot points based on quantiles of the original data)
  mc$formula <- mc_stored$formula
  mc$formula <- use_these_vars(y_mod_stuff$mod, 
                               c(y_mod_stuff$vars$formvars$fixed[-1], y_mod_stuff$vars$formvars$random[-1]), 
                               c(y_mod_stuff$vars$predvars$fixed[-1], y_mod_stuff$vars$predvars$random[-1]))
  mc$control <- mc_stored$control # revert to original control args    
  
  # Identify row in longitudinal data closest to event time or quadrature point
  #   NB if the quadrature point is earlier than the first observation time, 
  #   then covariates values are carried back to avoid missing values.
  #   In any other case, the observed covariates values from the most recent 
  #   observation time preceeding the quadrature point are carried forward to 
  #   represent the covariate value(s) at the quadrature point. (To avoid 
  #   missingness there is no limit on how far forwards or how far backwards 
  #   covariate values can be carried). If no time varying covariates are 
  #   present in the longitudinal submodel (other than the time variable) 
  #   then nothing is carried forward or backward.    
  mfq <- rolling_merge(data = mf, ids = e_flist, times = quadtimes)
  mod <- handle_glFormula(mc = mc, data = mfq, y_mod_stuff = y_mod_stuff)
  xtemp_eta <- mod$xtemp
  group_eta <- mod$group
  
  # If association structure is based on slope, then calculate design 
  # matrices under a time shift of epsilon
  sel_slope <- grep("etaslope|muslope", rownames(assoc))
  if (any(unlist(assoc[sel_slope,]))) {
    mfq_eps <- mfq
    mfq_eps[[time_var]] <- mfq_eps[[time_var]] + eps
    mod_eps <- handle_glFormula(mc = mc, data = mfq_eps, y_mod_stuff = y_mod_stuff)
    xtemp_eps <- mod_eps$xtemp
    group_eps <- mod_eps$group
  } else xtemp_eps <- group_eps <- NULL 
  
  # If association structure is based on a time lag, then calculate design 
  # matrices under the specified time lag
  sel_lag <- grep("etalag|mulag", rownames(assoc))
  if (any(unlist(assoc[sel_lag,]))) {
    quadtimes_lag <- lapply(quadtimes, function(x, lag) {
      newtimes <- x - lag
      newtimes[newtimes < 0] <- 0.0  # use baseline where lagged t is before baseline
      newtimes
    }, lag = ifelse(length(assoc["which_lag",][[m]]), assoc["which_lag",][[m]], 0))
    mfq_lag <- rolling_merge(data = mf, ids = e_flist, times = quadtimes_lag)
    mod_lag <- handle_glFormula(mc = mc, data = mfq_lag, y_mod_stuff = y_mod_stuff)
    xtemp_lag <- mod_lag$xtemp
    group_lag <- mod_lag$group
  } else xtemp_lag <- group_lag <- NULL
  
  # If association structure is based on area under the marker trajectory, then 
  # calculate design matrices at the subquadrature points
  sel_auc <- grep("etaauc|muauc", rownames(assoc))
  if (any(unlist(assoc[sel_auc,]))) {
    # Return a design matrix that is (quadnodes * auc_quadnodes * Npat) rows
    quadtimes_auc <- lapply(quadtimes, function(x) 
      unlist(lapply(x, function(y) 
        lapply(auc_quadpoints$points, unstandardise_quadpoints, 0, y))))
    ids2 <- rep(e_flist, each = auc_quadnodes)
    mfq_auc <- rolling_merge(data = mf, ids = ids2, times = quadtimes_auc)
    mod_auc <- handle_glFormula(mc = mc, data = mfq_auc, y_mod_stuff = y_mod_stuff)
    xtemp_auc <- mod_auc$xtemp
    group_auc <- mod_auc$group
  } else xtemp_auc <- group_auc <- NULL
  
  # If association structure is based on interactions with data, then calculate 
  # the design matrix which will be multiplied by etavalue, etaslope, muvalue or muslope
  sel_data <- grep("_data", rownames(assoc), value = TRUE)
  xq_data <- sapply(sel_data,
    function(x) {
      fm <- assoc["which_formulas",][[m]][[x]]
      if (length(fm)) {
        vars <- rownames(attr(terms.formula(fm), "factors"))
        if (is.null(vars))
          stop(paste0("No variables found in the formula specified for the '", x,
                      "' association structure.", call. = FALSE))
        ff <- ~ foo + bar
        gg <- parse(text = paste("~", paste(c(id_var, time_var), collapse = "+")))[[1L]]
        ff[[2L]][[2L]] <- fm[[2L]]
        ff[[2L]][[3L]] <- gg[[2L]]
        oldcall <- getCall(y_mod_stuff$mod)
        naa <- oldcall$na.action
        subset <- if (is.null(dataAssoc)) eval(oldcall$subset) else NULL               
        df <- if (is.null(dataAssoc)) eval(oldcall$data) else dataAssoc
        sel <- which(!vars %in% colnames(df))
        if (length(sel))
          stop(paste0("The following variables were specified in the formula for the '", x,
                      "' association structure, but they cannot be found in the data: ", 
                      paste0(vars, collapse = ", ")))
        mf2 <- eval(call("model.frame", ff, data = df, subset = subset, 
                         na.action = naa), envir = environment(y_mod_stuff$mod))
        mf2 <- data.table::data.table(mf2, key = c(id_var, time_var))
        mf2[[time_var]] <- as.numeric(mf2[[time_var]])
        mf2q <- get_quadrature_frame(data = mf2, ids = e_flist, times = quadtimes)
        xq <- stats::model.matrix(fm, data = mf2q)
        if ("(Intercept)" %in% colnames(xq)) xq <- xq[, -1L, drop = FALSE]
        if (!ncol(xq))
          stop(paste0("Bug found: A formula was specified for the '", x, "' association ", 
                      "structure, but the resulting design matrix has no columns."), call. = FALSE)
      } else xq <- matrix(0, length(unlist(quadtimes)), 0)
      xq
    }, simplify = FALSE, USE.NAMES = TRUE)
  K_data <- sapply(xq_data, ncol)
  xtemp_data <- do.call(cbind, xq_data)
  
  nlist(xtemp_eta, group_eta, xtemp_eps, group_eps, xtemp_lag, group_lag,
        xtemp_auc, group_auc, xq_data, xtemp_data, K_data, fr)  
}

# Carry out a rolling merge
#
# @param data A data.table with a set key corresponding to ids and times
# @param ids A vector of ids to merge against
# @param times A vector of (new) times to merge against
# @return A data.table formed by a merge of ids, times, and the closest 
#   preceding (in terms of times) rows in data
rolling_merge <- function(data, ids, times) {
  do.call(rbind, lapply(times, FUN = function(x) 
    data[data.table::SJ(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)]))  
}

# Evaluate a glFormula call and return model components
# 
# @param mc A glFormula call
# @param data A data frame to substitute into the data argument of the call
# @param y_mod_stuff A named list, returned by a call to handle_glmod (and
#   containing an indicator of whether the original longitudinal submodel had
#   an intercept term)
handle_glFormula <- function(mc, data, y_mod_stuff) { 
  mc$data <- data
  mod    <- eval(mc, parent.frame())
  x      <- as.matrix(mod$X)
  xtemp  <- if (y_mod_stuff$has_intercept) x[, -1L, drop = FALSE] else x  
  xtemp  <- sweep(xtemp, 2, y_mod_stuff$xbar, FUN = "-")
  group  <- mod$reTrms    
  nlist(xtemp, group)
}   

# Function to calculate the number of association parameters in the model
#
# @param assoc A list of length M with information about the association structure
#   type for each submodel, returned by an mapply call to validate_assoc
# @param a_mod_stuff A list of length M with the design matrices related to
#   the longitudinal submodels in the GK quadrature, returned by an mapply 
#   call to handle_a_mod
# @return Integer indicating the number of association parameters in the model 
get_num_assoc_pars <- function(assoc, a_mod_stuff) {
  sel1 <- c("etavalue", "etaslope", "etalag", "etaauc", 
            "muvalue", "muslope", "mulag", "muauc")
  sel2 <- c("which_b_zindex", "which_coef_zindex")
  sel3 <- c("which_interactions")
  K1 <- sum(as.integer(assoc[sel1,]))
  K2 <- length(unlist(assoc[sel2,]))
  K3 <- length(unlist(assoc[sel3,]))
  K4 <- sum(fetch_(a_mod_stuff, "K_data"))
  K1 + K2 + K3 + K4
}


#--------------- Functions related to prior weights

# Check the prior weights argument
#
# @param weights The data frame passed via the weights argument
# @param id_var The name of the ID variable
check_arg_weights <- function(weights, id_var) {
  
  # Check weights are an appropriate data frame
  if ((!is.data.frame(weights)) || (!ncol(weights) == 2))
    stop("'weights' argument should be a data frame with two columns: the first ",
         "containing patient IDs, the second containing their corresponding ",
         "weights.", call. = FALSE)
  if (!id_var %in% colnames(weights))
    stop("The data frame supplied in the 'weights' argument should have a ",
         "column named ", id_var, call. = FALSE)
  weight_var <- setdiff(colnames(weights), id_var)
  
  # Check weights are positive and numeric
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
}

# Return the vector of prior weights for one of the submodels
#
# @param mod_stuff A named list with elements: y, flist, ord
# @param weights The data frame passed via the weights argument
# @param id_var The name of the ID variable
handle_weights <- function(mod_stuff, weights, id_var) {
  
  is_glmod <- (is.null(mod_stuff$eventtime))
  
  # No weights provided by user
  if (is.null(weights)) {
    len <- if (is_glmod) length(mod_stuff$y) else length(mod_stuff$eventtime)
    return(rep(0.0, len)) 
  }

  # Check for IDs with no weight supplied
  weights[[id_var]] <- factor(weights[[id_var]])
  ids <- if (is_glmod) mod_stuff$flist[[id_var]] else factor(mod_stuff$flist)
  sel <- which(!ids %in% weights[[id_var]])
  if (length(sel)) {
    if (length(sel) > 30L) sel <- sel[1:30]
    stop(paste0("The following patient IDs are used in fitting the model, but ",
                "do not have weights supplied via the 'weights' argument: ",
                paste(ids[sel], collapse = ", ")), call. = FALSE)
  }
  
  # Obtain length and ordering of weights vector using flist
  wts_df  <- merge(data.frame(id = ids), weights, by.x = "id", by.y = id_var, sort = FALSE)
  wts_var <- setdiff(colnames(weights), id_var)
  wts     <- wts_df[[wts_var]]
  
  # Reorder weights if bernoulli
  ord <- mod_stuff[["ord"]]
  if (!is.null(ord)) wts <- wts[ord]

  wts
}


#--------------- Functions related to priors

# Autoscaling of priors
#
# @param prior_stuff A named list returned by a call to handle_glm_prior
# @param mod_stuff A named list returned by a call to either handle_glmod or
#   handle_coxmod
# @param QR A logical specifying whether QR decomposition is used for the x matrix
# @param use_x A logical specifying whether to autoscale the priors based on
#   the standard deviations of the predictor variables
# @param min_prior_scale The minimum allowed for prior scales
# @return A named list of the same structure as returned by handle_glm_prior
autoscale_prior <- function(prior_stuff, mod_stuff, QR, use_x = FALSE, min_prior_scale = 1e-12) {
  
  is_glmod <- is.null(mod_stuff$eventtime)
  
  if (is_glmod && mod_stuff$is_gaussian) {
    ss <- sd(mod_stuff$y)
    if (prior_stuff$prior_dist > 0L && prior_stuff$prior_autoscale)
      prior_stuff$prior_scale <- ss * prior_stuff$prior_scale
  }
  
  if (use_x) {
    if (!QR && prior_stuff$prior_dist > 0L && prior_stuff$prior_autoscale) {
      prior_stuff$prior_scale <- pmax(min_prior_scale, prior_stuff$prior_scale / 
                                        apply(mod_stuff$xtemp, 2L, FUN = function(x) {
                                          num.categories <- length(unique(x))
                                          x.scale <- 1
                                          if (num.categories == 2) {
                                            x.scale <- diff(range(x))
                                          } else if (num.categories > 2) {
                                            x.scale <- sd(x)
                                          }
                                          return(x.scale)
                                        }))
    }      
  }
  prior_stuff$prior_scale <- 
    as.array(pmin(.Machine$double.xmax, prior_stuff$prior_scale))
  
  prior_stuff
}

# Get the required number of (local) horseshoe parameters for a specified prior type
#
# @param prior_dist An integer indicating the type of prior distribution: 
#   where 1L == normal, 2L == t, 3L == hs, 4L == hs_plus
get_nvars_for_hs <- function(prior_dist) {
  if      (prior_dist <= 2L) return(0L) 
  else if (prior_dist == 3L) return(2L) 
  else if (prior_dist == 4L) return(4L)
  else return(0L)
}


#--------------- Functions related to generating initial values

# Create a function that can be used to generate the model-based initial values for Stan
#
# @param y_mod_stuff A list, with each element containing the list object returned by
#   a call to the handle_glmod function
# @param e_mod_stuff A list object returned by a call to the handle_coxmod function
# @param standata The data list that will be passed to Stan
generate_init_function <- function(y_mod_stuff, e_mod_stuff, standata) {
  
  # Initial values for intercepts, coefficients and aux parameters
  y_est <- lapply(y_mod_stuff, function(x) summary(x$mod)$coefficients[, "Estimate"])
  y_xbar          <- fetch(y_mod_stuff, "xbar")
  y_gamma         <- lapply(seq_along(y_mod_stuff), function(m) 
    return_intercept(y_est[[m]]) - y_xbar[[m]] %*% drop_intercept(y_est[[m]]))
  y_gamma_unbound <- y_gamma[as.logical(standata$y_has_intercept_unbound)]
  y_gamma_lobound <- y_gamma[as.logical(standata$y_has_intercept_lobound)]
  y_gamma_upbound <- y_gamma[as.logical(standata$y_has_intercept_upbound)]
  y_beta          <- lapply(y_est, drop_intercept)
  e_beta          <- e_mod_stuff$mod$coef
  y_aux           <- lapply(y_mod_stuff, function(x) sigma(x$mod))[as.logical(standata$y_has_aux)]
  e_aux           <- if (standata$basehaz_type == 1L) runif(1, 0.5, 3) else rep(0, standata$basehaz_df)
  y_z_beta        <- standardise_coef(unlist(y_beta), standata$y_prior_mean,         standata$y_prior_scale)
  e_z_beta        <- standardise_coef(e_beta,         standata$e_prior_mean,         standata$e_prior_scale) 
  y_aux_unscaled  <- standardise_coef(unlist(y_aux),  standata$y_prior_mean_for_aux, standata$y_prior_scale_for_aux)
  e_aux_unscaled  <- standardise_coef(e_aux,          standata$e_prior_mean_for_aux, standata$e_prior_scale_for_aux)
  b_Cov           <- lapply(y_mod_stuff, function(x) lme4::VarCorr(x$mod)[[1L]])
  sel             <- sapply(y_mod_stuff, function(x) length(lme4::VarCorr(x$mod)) > 1L)
  if (any(sel)) stop("Model-based initial values cannot yet be used with more ",
                     "than one clustering level.", call. = FALSE)
 
  # Initial values for random effects distribution
  scale <- standata$scale
  t     <- standata$t
  p     <- standata$p

  # Cholesky decomp of b_Cov combined across all submodel
  L_b_Cov         <- t(suppressWarnings(chol(as.matrix(Matrix::bdiag(b_Cov)), pivot = TRUE))) 
  diag_L_b_Cov    <- diag(L_b_Cov)
  
  # Dimensions
  len_z_T <- 0
  for (i in 1:t) {
    if (p[i] > 2) 
      for (j in 3:p[i]) len_z_T <- len_z_T + p[i] - 1;
  }
  len_rho <- ifelse((sum(p) - t) > 0, sum(p) - t, 0)

  # Construct initial values for theta_L matrix
  # ** Much room for improvement here! **
  tau              <- c()
  rho              <- c()
  normalised_zetas <- c()
  rho_mark         <- 1
  for (i in 1:t) {
    trace <- sum(diag_L_b_Cov)  # equal to variance of RE if only one RE
    normalised_zetas_tmp <- diag_L_b_Cov / trace  # equal to 1 if only one random effect
    tau[i] <- (sqrt(trace / p[i])) / scale[i]
    std_dev1 <- sqrt(L_b_Cov[1,1])
    if (p[i] > 1) {
      std_dev2 <- sqrt(L_b_Cov[2,2])
      T21 <- L_b_Cov[2,1] / std_dev2
      rho[rho_mark] <- (T21 + 1) / 2
      rho_mark <- rho_mark + 1
      normalised_zetas <- c(normalised_zetas, normalised_zetas_tmp)
    }
  }
  
  # Function to generate model based initial values
  model_based_inits <- Filter(function(x) (!is.null(x)), list(
    y_gamma_unbound = array_else_double(y_gamma_unbound),
    y_gamma_lobound = array_else_double(y_gamma_lobound),
    y_gamma_upbound = array_else_double(y_gamma_upbound),
    y_z_beta        = array_else_double(y_z_beta),
    y_aux_unscaled  = array_else_double(y_aux_unscaled),
    e_z_beta        = array_else_double(e_z_beta),
    e_aux_unscaled  = array_else_double(e_aux_unscaled),
    e_gamma         = array_else_double(rep(0, standata$e_has_intercept)),
    a_z_beta        = array_else_double(rep(0, standata$a_K)),
    z_b             = array_else_double(runif(standata$len_b, -0.5, 0.5)),
    y_global        = array_else_double(runif(get_nvars_for_hs(standata$y_prior_dist))),
    e_global        = array_else_double(runif(get_nvars_for_hs(standata$e_prior_dist))),
    a_global        = array_else_double(runif(get_nvars_for_hs(standata$a_prior_dist))),
    y_local         = matrix_of_uniforms(nrow = get_nvars_for_hs(standata$y_prior_dist), ncol = standata$sum_y_K),
    e_local         = matrix_of_uniforms(nrow = get_nvars_for_hs(standata$e_prior_dist), ncol = standata$e_K),
    a_local         = matrix_of_uniforms(nrow = get_nvars_for_hs(standata$a_prior_dist), ncol = standata$a_K),
    y_S             = if (standata$y_prior_dist %in% c(5,6)) matrix(rep(1, standata$sum_y_K), 1, standata$sum_y_K) else matrix(0, 0, standata$sum_y_K),
    e_S             = if (standata$e_prior_dist %in% c(5,6)) matrix(rep(1, standata$e_K), 1, standata$e_K) else matrix(0, 0, standata$e_K),
    a_S             = if (standata$a_prior_dist %in% c(5,6)) matrix(rep(1, standata$a_K), 1, standata$a_K) else matrix(0, 0, standata$a_K),
    y_one_over_lambda = if (standata$y_prior_dist == 6) as.array(1) else as.array(double(0)), 
    e_one_over_lambda = if (standata$e_prior_dist == 6) as.array(1) else as.array(double(0)), 
    a_one_over_lambda = if (standata$a_prior_dist == 6) as.array(1) else as.array(double(0)),
    z_T             = array_else_double(rep(sqrt(1 / len_z_T), len_z_T)),
    rho             = array_else_double(rep(1 / (len_rho + 1), len_rho)),
    zeta            = array_else_double(normalised_zetas),
    tau             = array_else_double(tau)))
  
  return(function() model_based_inits)
}


#--------------- Functions related to standata and sampling

# Set arguments for sampling for stan_jm
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# *Note that this differs from the set_sampling_args function in that
# it uses different default numbers of iterations and chains, and 
# different default adapt_delta and max_treedepth. Using the shorter 
# treedepth is to stop the sampler trailing off during early iterations. 
# This can drastically reduce the model estimation time, and in most
# examples using a shorter treedepth hasn't compromised the sampler
# at later interations (ie, at later iterations the sampler doesn't
# hit the maximum treedepth).
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
  if (!"refresh" %in% unms) args$refresh <- args$iter / 25
  if (!"save_warmup" %in% unms) args$save_warmup <- TRUE  
  
  return(args)
}  

#--------------- Miscellaneous and helper functions

# Check argument input type is ok, and return as a list
validate_arg <- function(arg, type, null_ok = FALSE, 
                         validate_length = NULL, broadcast = TRUE) {
  if (is.null(arg)) { # input type NULL, check if ok and return list or error
    if (null_ok) arg <- list(arg) else STOP_arg(arg, type, null_ok = null_ok)
  } else if (any(sapply(type, function(x) is(arg, x)))) { # input type ok, return as list
    arg <- list(arg)
  } else if (is(arg, "list")) { # input list, check each element
    check <- sapply(arg, function(x) if (null_ok) 
      (is.null(x) || any(sapply(type, function(y) is(x, y)))) else 
        (any(sapply(type, function(y) is(x, y)))))
    if (!all(check)) STOP_arg(arg, type, null_ok = null_ok)
  } else {
    STOP_arg(arg, type, null_ok = null_ok)
  }
  if (!is.null(validate_length)) {
    if (length(arg) == 1L && validate_length > 1) {
      arg <- if (broadcast) rep(arg, times = validate_length) else obj
    } else if (!length(arg) == validate_length) {
      stop(paste0(arg, " is a list of the incorrect length."))
    }
  }
  arg
}

# Check if the user input a list of priors for the longitudinal
# submodel, and if not, then return the appropriate list
#
# @param priorarg The user input to the priorLong* argument in the stan_jm call
# @param M An integer specifying the number of longitudinal submodels
maybe_broadcast_priorarg <- function(priorarg, M) {
  if (is.null(priorarg)) {
    return(rep(list(NULL), M))
  } else if ("dist" %in% names(priorarg)) {
    return(rep(list(priorarg), M))
  } else if (is.list(priorarg) && (length(priorarg) == M)) {
    return(priorarg)
  } else {
    nm <- deparse(substitute(priorarg))
    stop(nm, "appears to provide prior information separately for the different ",
         "longitudinal submodels, but the list is of the incorrect length.")
  }
}

# Error message when the argument contains an object of the incorrect type
STOP_arg <- function(arg, type, null_ok = FALSE) {
  stop(paste0(deparse(substitute(arg)), " should be ", ifelse(null_ok, "NULL, ", ""),
              "a ", type, " object or a list of ", type, " objects.")) 
}

# Return error msg if both elements of the object are TRUE
STOP_combination_not_allowed <- function(object, x, y) {
  if (object[[x]] && object[[y]])
    stop("In ", deparse(substitute(object)), ", '", x, "' and '", y,
         "' cannot be specified together", call. = FALSE)
}

# Return a list (or vector if unlist = TRUE) which
# contains the embedded elements in list x named y 
fetch <- function(x, y, unlist = FALSE) {
  ret <- lapply(x, `[[`, y)
  if (unlist) unlist(ret) else ret
}
# Wrapper for using fetch with unlist = TRUE
fetch_ <- function(x, y) {
  fetch(x, y, unlist = TRUE)
}
# Wrapper for using fetch with unlist = TRUE 
# and returning array
fetch_array <- function(x, y) {
  as.array(fetch(x, y, unlist = TRUE))
}

# Drop intercept from a vector of named coefficients
drop_intercept <- function(x) { 
  sel <- which("(Intercept)" %in% names(x))
  if (length(sel)) x[-sel] else x
}

# Return intercept from a vector of named coefficients
return_intercept <- function(x) {
  sel <- which("(Intercept)" %in% names(x))
  if (length(sel)) x[sel] else NULL
}

# Standardise a coefficient
standardise_coef <- function(x, location = 0, scale = 1)
  (x - location) / scale

# Return a one-dimensional array or an empty numeric
array_else_double <- function(x)
  if (!length(x)) double(0) else as.array(unlist(x))

# Return a matrix of uniform random variables or an empty matrix
matrix_of_uniforms <- function(nrow = 0, ncol = 0)
  if (nrow == 0 || ncol == 0) matrix(0,0,0) else matrix(runif(nrow * ncol), nrow, ncol)




