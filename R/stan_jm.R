# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

#' Bayesian joint longitudinal and time-to-event models via Stan
#' 
#' Fits a shared parameter joint model for longitudinal and time-to-event 
#' (e.g. survival) data under a Bayesian framework using Stan.
#' 
#' @export
#' @templateVar pkg stats
#' @templateVar pkgfun glm
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
#'   be either a single data frame which contains the data for all 
#'   longitudinal submodels, or it can be a list of data frames where each
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
#'   be used include: "etavalue" (the default); "etaslope"; "etaauc"; 
#'   "muvalue"; "muslope"; "muauc"; "shared_b"; "shared_coef"; or "null". 
#'   These are described in the \strong{Details} section below. For a multivariate 
#'   joint model, different association structures can optionally be used for 
#'   each longitudinal submodel by specifying a list of character
#'   vectors, with each element of the list specifying the desired association 
#'   structure for one of the longitudinal submodels. Specifying \code{assoc = NULL}
#'   will fit a joint model with no association structure (equivalent  
#'   to fitting separate longitudinal and time-to-event models). It is also 
#'   possible to include interaction terms between the association term 
#'   ("etavalue", "etaslope", "muvalue", "muslope") and observed data/covariates. 
#'   It is also possible, when fitting a multivariate joint model, to include 
#'   interaction terms between the association terms ("etavalue" or "muvalue") 
#'   corresponding to the different longitudinal outcomes. See the 
#'   \strong{Details} section as well as the \strong{Examples} below.
#' @param lag_assoc A non-negative scalar specifying the time lag that should be
#'   used for the association structure. That is, the hazard of the event at 
#'   time \emph{t} will be assumed to be associated with the value/slope/auc of 
#'   the longitudinal marker at time \emph{t-u}, where \emph{u} is the time lag.
#'   If fitting a multivariate joint model, then a different time lag can be used
#'   for each longitudinal marker by providing a numeric vector of lags, otherwise
#'   if a scalar is provided then the specified time lag will be used for all 
#'   longitudinal markers. Note however that only one time lag  can be specified 
#'   for linking each longitudinal marker to the 
#'   event, and that that time lag will be used for all association structure
#'   types (e.g. \code{"etavalue"}, \code{"etaslope"}, \code{"etaauc"}, 
#'   \code{"muvalue"}, etc) that are specified for that longitudinal marker in
#'   the \code{assoc} argument.
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
#' @param weights Experimental and should be used with caution. The 
#'   user can optionally supply a 2-column data frame containing a set of
#'   'prior weights' to be used in the estimation process. The data frame should
#'   contain two columns: the first containing the IDs for each individual, and 
#'   the second containing the corresponding weights. The data frame should only
#'   have one row for each individual; that is, weights should be constant 
#'   within individuals.
#' @param init The method for generating the initial values for the MCMC.
#'   The default is \code{"prefit"}, which uses those obtained from 
#'   fitting separate longitudinal and time-to-event models prior  
#'   to fitting the joint model. The separate models are estimated using the
#'   \code{\link[lme4]{glmer}} and \code{\link[survival]{coxph}} functions.
#'   Parameters that cannot be obtained from 
#'   fitting separate longitudinal and time-to-event models are initialised 
#'   at 0. This often provides reasonable initial values which should aid the 
#'   MCMC sampler, but may not work for some complex or heavily constrained 
#'   models such as those with multiple longitudinal submodels or multilevel
#'   clustering. An alternative is to \code{"prefit_vb"} which will obtain 
#'   initial values for the longitudinal submodels by estimating
#'   a multivariate GLM using the \code{\link{stan_mvmer}} function with 
#'   \code{algorithm = "meanfield"}. This may help with initial values for some
#'   complex models. Note that it is recommended that any final analysis should 
#'   be performed with several MCMC chains each initiated from a different
#'   set of initial values; this can be obtained by setting
#'   \code{init = "random"}. Other possibilities for specifying \code{init}
#'   are those described for \code{\link[rstan]{stan}}.  
#' @param priorLong,priorEvent,priorAssoc The prior distributions for the 
#'   regression coefficients in the longitudinal submodel(s), event submodel,
#'   and the association parameter(s). Can be a call to one of the various functions 
#'   provided by \pkg{rstanarm} for specifying priors. The subset of these functions 
#'   that can be used for the prior on the coefficients can be grouped into several 
#'   "families":
#'   
#'   \tabular{ll}{
#'     \strong{Family} \tab \strong{Functions} \cr 
#'     \emph{Student t family} \tab \code{normal}, \code{student_t}, \code{cauchy} \cr 
#'     \emph{Hierarchical shrinkage family} \tab \code{hs}, \code{hs_plus} \cr 
#'     \emph{Laplace family} \tab \code{laplace}, \code{lasso} \cr
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
#' @param max_treedepth A positive integer specifying the maximum treedepth 
#'   for the non-U-turn sampler. See the \code{control} argument in 
#'   \code{\link[rstan]{stan}}.
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
#'   \code{\link[lme4]{glmer}}. If a multivariate joint model is specified (by
#'   providing a list of formulas in the \code{formulaLong} argument), then
#'   the multivariate longitudinal submodel consists of a multivariate generalized  
#'   linear model (GLM) with group-specific terms that are assumed to be correlated
#'   across the different GLM submodels. That is, within
#'   a grouping factor (for example, patient ID) the group-specific terms are
#'   assumed to be correlated across the different GLM submodels. It is 
#'   possible to specify a different outcome type (for example a different
#'   family and/or link function) for each of the GLM submodels, by providing
#'   a list of \code{\link[stats]{family}} objects in the \code{family} 
#'   argument. \cr
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
#'   Time-varying covariates are allowed in both the 
#'   longitudinal and event submodels. These should be specified in the data 
#'   in the same way as they normally would when fitting a separate 
#'   longitudinal model using \code{\link[lme4]{lmer}} or a separate 
#'   time-to-event model using \code{\link[survival]{coxph}}. These time-varying
#'   covariates should be exogenous in nature, otherwise they would perhaps 
#'   be better specified as an additional outcome (i.e. by including them as an 
#'   additional longitudinal outcome in the joint model). \cr
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
#'   \subsection{Association structures}{
#'   The association structure for the joint model can be based on any of the 
#'   following parameterisations: 
#'     \itemize{
#'       \item current value of the linear predictor in the 
#'         longitudinal submodel (\code{"etavalue"}) 
#'       \item first derivative (slope) of the linear predictor in the 
#'         longitudinal submodel (\code{"etaslope"}) 
#'       \item the area under the curve of the linear predictor in the 
#'         longitudinal submodel (\code{"etaauc"}) 
#'       \item current expected value of the longitudinal submodel 
#'         (\code{"muvalue"})
#'       \item the area under the curve of the expected value from the 
#'         longitudinal submodel (\code{"muauc"})
#'       \item shared individual-level random effects (\code{"shared_b"}) 
#'       \item shared individual-level random effects which also incorporate 
#'         the corresponding fixed effect as well as any corresponding 
#'         random effects for clustering levels higher than the individual)
#'         (\code{"shared_coef"})
#'       \item interactions between association terms and observed data/covariates
#'         (\code{"etavalue_data"}, \code{"etaslope_data"}, \code{"muvalue_data"}, 
#'         \code{"muslope_data"}). These are described further below.
#'       \item interactions between association terms corresponding to different 
#'         longitudinal outcomes in a multivariate joint model 
#'         (\code{"etavalue_etavalue(#)"}, \code{"etavalue_muvalue(#)"},
#'         \code{"muvalue_etavalue(#)"}, \code{"muvalue_muvalue(#)"}). These
#'         are described further below.      
#'       \item no association structure (equivalent to fitting separate 
#'         longitudinal and event models) (\code{"null"} or \code{NULL}) 
#'     }
#'   More than one association structure can be specified, however,
#'   not all possible combinations are allowed.   
#'   Note that for the lagged association structures baseline values (time = 0) 
#'   are used for the instances 
#'   where the time lag results in a time prior to baseline. When using the 
#'   \code{"etaauc"} or \code{"muauc"} association structures, the area under
#'   the curve is evaluated using Gauss-Kronrod quadrature with 15 quadrature 
#'   nodes. By default, \code{"shared_b"} and \code{"shared_coef"} contribute 
#'   all random effects to the association structure; however, a subset of the 
#'   random effects can be chosen by specifying their indices between parentheses 
#'   as a suffix, for example, \code{"shared_b(1)"} or \code{"shared_b(1:3)"} or 
#'   \code{"shared_b(1,2,4)"}, and so on. \cr
#'   \cr 
#'   In addition, several association terms (\code{"etavalue"}, \code{"etaslope"},
#'   \code{"muvalue"}, \code{"muslope"}) can be interacted with observed 
#'   data/covariates. To do this, use the association term's main handle plus a
#'   suffix of \code{"_data"} then followed by the model matrix formula in 
#'   parentheses. For example if we had a variable in our dataset for gender 
#'   named \code{sex} then we might want to obtain different estimates for the 
#'   association between the current slope of the marker and the risk of the 
#'   event for each gender. To do this we would specify 
#'   \code{assoc = "etaslope_data(~ sex)"}. \cr
#'   \cr
#'   It is also possible, when fitting  a multivariate joint model, to include 
#'   interaction terms between the association terms themselves (this only
#'   applies for interacting \code{"etavalue"} or \code{"muvalue"}). For example, 
#'   if we had a joint model with two longitudinal markers, we could specify 
#'   \code{assoc = list(c("etavalue", "etavalue_etavalue(2)"), "etavalue")}.
#'   The first element of list says we want to use the value of the linear
#'   predictor for the first marker, as well as it's interaction with the
#'   value of the linear predictor for the second marker. The second element of 
#'   the list says we want to also include the expected value of the second marker 
#'   (i.e. as a "main effect"). Therefore, the linear predictor for the event 
#'   submodel would include the "main effects" for each marker as well as their
#'   interaction. \cr
#'   \cr
#'   There are additional examples in the \strong{Examples} section below.
#'   }
#' 
#' @return A \link[=stanmvreg-objects]{stanmvreg} object is returned.
#' 
#' @seealso \code{\link{stanmvreg-objects}}, \code{\link{stanmvreg-methods}}, 
#'   \code{\link{print.stanmvreg}}, \code{\link{summary.stanmvreg}},
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
#'               dataLong = pbcLong,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv,
#'               time_var = "year",
#'               # this next line is only to keep the example small in size!
#'               chains = 1, cores = 1, seed = 12345, iter = 1000)
#' summary(f1) 
#'         
#' #####
#' # Univariate joint model, with association structure based on the 
#' # current value and slope of the linear predictor
#' f2 <- stan_jm(formulaLong = logBili ~ year + (year | id), 
#'               dataLong = pbcLong,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv,
#'               assoc = c("etavalue", "etaslope"),
#'               time_var = "year",
#'               chains = 1, cores = 1, seed = 12345, iter = 1000)
#' summary(f2)  
#' 
#' #####
#' # Univariate joint model, with association structure based on the 
#' # lagged value of the linear predictor, where the lag is 2 time 
#' # units (i.e. 2 years in this example)
#' f3 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv,
#'               time_var = "year",
#'               assoc = "etavalue", lag_assoc = 2,
#'               chains = 1, cores = 1, seed = 12345, iter = 1000)
#' summary(f3) 
#' 
#' #####
#' # Univariate joint model, where the association structure includes 
#' # interactions with observed data. Here we specify that we want to use 
#' # an association structure based on the current value of the linear 
#' # predictor from the longitudinal submodel (i.e. "etavalue"), but we 
#' # also want to interact this with the treatment covariate (trt) from
#' # pbcLong data frame, so that we can estimate a different association 
#' # parameter (i.e. estimated effect of log serum bilirubin on the log 
#' # hazard of death) for each treatment group
#' f4 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv,
#'               time_var = "year", chains = 1,
#'               assoc = c("etavalue", "etavalue_data(~ trt)"))
#' 
#' ######
#' # Multivariate joint model, with association structure based 
#' # on the current value of the linear predictor in the first longitudinal 
#' # submodel and shared random intercept from the second longitudinal 
#' # submodel only (which is the first random effect in that submodel
#' # and is therefore indexed using the '(1)' suffix in the code below)
#' mv1 <- stan_jm(
#'         formulaLong = list(
#'           logBili ~ year + (1 | id), 
#'           albumin ~ sex + year + (year | id)),
#'         dataLong = pbcLong,
#'         formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'         dataEvent = pbcSurv,
#'         assoc = list("etavalue", "shared_b(1)"), 
#'         time_var = "year",
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' summary(mv1)
#' 
#' # To include both the random intercept and random slope in the shared 
#' # random effects association structure for the second longitudinal 
#' # submodel, we could specify the following:
#' #   update(mv1, assoc = list("etavalue", "shared_b")
#' # which would be equivalent to:  
#' #   update(mv1, assoc = list("etavalue", "shared_b(1,2)")
#' # or:
#' #   update(mv1, assoc = list("etavalue", "shared_b(1:2)")     
#' 
#' #####
#' # Multivariate joint model, where the association structure is formed by 
#' # including the expected value of each longitudinal marker (logBili and 
#' # albumin) in the linear predictor of the event submodel, as well as their 
#' # interaction effect (i.e. the interaction between the two "etavalue" terms). 
#' # Note that whether such an association structure based on a marker by 
#' # marker interaction term makes sense will depend on the context of your 
#' # application -- here we just show it for demostration purposes).
#' mv2 <- stan_jm(
#'         formulaLong = list(
#'           logBili ~ year + (1 | id), 
#'           albumin ~ sex + year + (year | id)),
#'         dataLong = pbcLong,
#'         formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'         dataEvent = pbcSurv,
#'         assoc = list(c("etavalue", "etavalue_etavalue(2)"), "etavalue"),
#'         time_var = "year", 
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' }
#'  
#' @import data.table
#' @importFrom lme4 lmerControl glmerControl glmer
#' 
stan_jm <- function(formulaLong, dataLong, formulaEvent, dataEvent, time_var, 
                    id_var, family = gaussian, assoc = "etavalue", 
                    lag_assoc = 0, dataAssoc,
                    basehaz = c("weibull", "bs", "piecewise"), basehaz_ops, 
                    quadnodes = 15, init = "prefit", 
                    na.action = getOption("na.action", "na.omit"), weights, 
                    offset, contrasts, ...,				          
                    priorLong = normal(), priorLong_intercept = normal(), 
                    priorLong_aux = cauchy(0, 5), priorEvent = normal(), 
                    priorEvent_intercept = normal(), priorEvent_aux = cauchy(0, 50),
                    priorAssoc = normal(), prior_covariance = decov(), prior_PD = FALSE, 
                    algorithm = c("sampling", "meanfield", "fullrank"), 
                    adapt_delta = NULL, max_treedepth = NULL, QR = FALSE, 
                    sparse = FALSE) {
  
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  # Set seed if specified
  dots <- list(...)
  if ("seed" %in% names(dots))
    set.seed(dots$seed)
  
  # Check for arguments not yet implemented
  if (!missing(dataAssoc))
    stop("'dataAssoc argument not yet implemented.")
  if (!missing(offset)) 
    stop("Offsets are not yet implemented for stan_jm")
  if (QR)               
    stop("QR decomposition not yet implemented for stan_jm")
  if (sparse)
    stop("'sparse' option is not yet implemented for stan_jm")
  if (missing(offset))      offset      <- NULL 
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  if (missing(weights))     weights     <- NULL
  if (missing(id_var))      id_var      <- NULL
  if (missing(dataAssoc))   dataAssoc   <- NULL 
  
  # Validate arguments
  basehaz   <- match.arg(basehaz)
  algorithm <- match.arg(algorithm)
  formulaLong <- validate_arg(formulaLong, "formula")
  M           <- length(formulaLong)
  dataLong    <- validate_arg(dataLong, "data.frame", null_ok = TRUE, 
                              validate_length = M, broadcast = TRUE)
  assoc       <- validate_arg(assoc, "character", null_ok = TRUE, 
                              validate_length = M, broadcast = TRUE)
  dataLong <- lapply(dataLong, as.data.frame)
  dataEvent <- as.data.frame(dataEvent)

  # Check family and link
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  if (!is(family, "list")) {
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
  calling_env <- parent.frame()
  call <- match.call(expand.dots = TRUE)    
  mc   <- match.call(expand.dots = FALSE)
  mc$time_var <- mc$id_var <- mc$assoc <- mc$lag_assoc <- 
    mc$dataAssoc <- mc$basehaz <- mc$basehaz_ops <-
    mc$df <- mc$knots <- mc$quadnodes <- NULL
  mc$priorLong <- mc$priorLong_intercept <- mc$priorLong_aux <-
    mc$priorEvent <- mc$priorEvent_intercept <- mc$priorEvent_aux <-
    mc$priorAssoc <- mc$prior_covariance <-  mc$prior_PD <- mc$algorithm <- 
    mc$scale <- mc$concentration <- mc$shape <- mc$init <- mc$adapt_delta <- 
    mc$max_treedepth <- mc$... <- mc$QR <- NULL
  mc$weights <- NULL 

  # Create call for longitudinal submodel  
  y_mc <- mc
  y_mc <- strip_nms(y_mc, "Long") 
  y_mc$formulaEvent <- y_mc$dataEvent <- NULL

  # Create call for each longitudinal submodel separately
  m_mc <- lapply(1:M, function(m, old_call, env) {
    new_call <- old_call
    fm     <- eval(old_call$formula, env)
    data   <- eval(old_call$data, env)
    family <- eval(old_call$family, env)
    new_call$formula <- if (is(fm, 'list')) fm[[m]] else fm
    new_call$data    <- if (is(data, 'list') && !inherits(data, 'data.frame')) data[[m]] else data
    new_call$family  <- if (is(family, 'list')) family[[m]] else family
    new_call
  }, old_call = y_mc, env = calling_env)

  # Create call for event submodel
  e_mc <- mc
  e_mc <- strip_nms(e_mc, "Event")
  e_mc$formulaLong <- e_mc$dataLong <- e_mc$family <- NULL
  
  # Is priorLong* already a list?
  priorLong           <- maybe_broadcast_priorarg(priorLong, M)
  priorLong_intercept <- maybe_broadcast_priorarg(priorLong_intercept, M)
  priorLong_aux       <- maybe_broadcast_priorarg(priorLong_aux, M)
    
  #--------------------------------
  # Data for longitudinal submodel
  #--------------------------------
  
  # Fit separate longitudinal submodels
  y_mod_stuff <- mapply(handle_glmod, m_mc, family, 
                        MoreArgs = list(supported_families = supported_families, 
                                        supported_links = supported_links, 
                                        sparse = sparse, env = calling_env), 
                        SIMPLIFY = FALSE)

  # Construct single cnms list for all longitudinal submodels
  cnms <- get_common_cnms(fetch(y_mod_stuff, "cnms"))
  cnms_nms <- names(cnms)
  
  # Additional error checks
  id_var <- check_id_var(id_var, fetch(y_mod_stuff, "cnms"), fetch(y_mod_stuff, "flist"))
  unique_id_list <- check_id_list(id_var, fetch(y_mod_stuff, "flist"))
    
  # Construct prior weights
  has_weights <- (!is.null(weights))
  if (has_weights) check_arg_weights(weights, id_var)
  y_weights <- lapply(y_mod_stuff, handle_weights, weights, id_var)
  
  #-------------------------
  # Data for event submodel
  #-------------------------
  
  if (!id_var %in% colnames(dataEvent))
    stop(paste0("Variable '", id_var, "' must be appear in dataEvent"), call. = FALSE)
  
  # Fit separate event submodel
  e_mod_stuff <- handle_coxmod(e_mc, quadnodes = quadnodes, id_var = id_var, 
                               unique_id_list = unique_id_list, sparse = sparse,
                               env = calling_env)
  
  # Construct prior weights
  e_weights <- handle_weights(e_mod_stuff, weights, id_var)

  # Baseline hazard
  ok_basehaz <- nlist("weibull", "bs", "piecewise")
  basehaz <- handle_basehaz(basehaz, basehaz_ops, ok_basehaz = ok_basehaz, 
                            eventtime = e_mod_stuff$eventtime, d = e_mod_stuff$d)
  
  # Incorporate intercept term if Weibull baseline hazard
  e_has_intercept <- e_mod_stuff$has_intercept <- (basehaz$type == 1L)

  #--------------------------------
  # Data for association structure
  #--------------------------------
  
  # Handle association structure
  # !!! NB if ordering is changed here, then must also change standata$has_assoc
  ok_assoc <- c("null", "etavalue","etaslope", "etaauc", "muvalue", 
                "muslope", "muauc", "shared_b", "shared_coef")
  ok_assoc_data         <- ok_assoc[c(2:3,5:6)]
  ok_assoc_interactions <- ok_assoc[c(2,5)]
  
  # Check lag
  if (length(lag_assoc) == 1L)
    lag_assoc <- rep(lag_assoc, M)
  if (!length(lag_assoc) == M)
    stop("'lag_assoc' should length 1 or length equal to the number of markers.")
  if (!is.numeric(lag_assoc))
    stop("'lag_assoc' must be numeric.")
  if (any(lag_assoc < 0))
    stop("'lag_assoc' must be non-negative.")
  
  assoc <- mapply(validate_assoc, assoc, y_mod_stuff, lag = lag_assoc,
                  MoreArgs = list(ok_assoc = ok_assoc, ok_assoc_data = ok_assoc_data,
                                  ok_assoc_interactions = ok_assoc_interactions, 
                                  id_var = id_var, M = M))
  assoc <- check_order_of_assoc_interactions(assoc, ok_assoc_interactions)
  colnames(assoc) <- paste0("Long", 1:M)

  # Return design matrices for evaluating longitudinal submodel quantities
  # at the quadrature points
  eps <- 1E-5 # time shift for numerically calculating deriv using one-sided diff
  auc_quadnodes <- 15L
  a_mod_stuff <- mapply(handle_assocmod, 1:M, m_mc, dataLong, y_mod_stuff, 
                        SIMPLIFY = FALSE, 
                        MoreArgs = list(id_list         = e_mod_stuff$flist, 
                                        times           = e_mod_stuff$quadtimes, 
                                        assoc           = assoc, 
                                        id_var          = id_var, 
                                        time_var        = time_var, 
                                        eps             = eps, 
                                        auc_quadnodes   = auc_quadnodes,
                                        dataAssoc       = dataAssoc,
                                        env             = calling_env))

  # Number of association parameters
  a_K <- get_num_assoc_pars(assoc, a_mod_stuff)
  
  #---------------------
  # Prior distributions
  #---------------------
  
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso")
  ok_intercept_dists <- ok_dists[1:3]
  ok_y_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  ok_e_aux_dists <- ok_dists[1:3]
  
  # Note: *_user_prior_*_stuff objects are stored unchanged for constructing 
  # prior_summary, while *_prior_*_stuff objects are autoscaled
  
  # Priors for longitudinal submodel(s)
  y_user_prior_stuff <- y_prior_stuff <- mapply(
    handle_glm_prior,
    priorLong,
    nvars = fetch(y_mod_stuff, "K"),
    link = fetch(family, "link"),
    MoreArgs = list(default_scale = 2.5, ok_dists = ok_dists), 
    SIMPLIFY = FALSE)

  y_user_prior_intercept_stuff <- y_prior_intercept_stuff <- mapply(
    handle_glm_prior,
    priorLong_intercept,
    link = fetch(family, "link"),
    MoreArgs = list(nvars = 1, default_scale = 10, ok_dists = ok_intercept_dists), 
    SIMPLIFY = FALSE)

  y_user_prior_aux_stuff <- y_prior_aux_stuff <- mapply(
    handle_glm_prior,
    priorLong_aux,
    MoreArgs = list(nvars = 1, default_scale = 5, link = NULL, ok_dists = ok_y_aux_dists), 
    SIMPLIFY = FALSE)  
  
  # Priors for event submodel
  e_user_prior_stuff <- e_prior_stuff <- 
    handle_glm_prior(
      priorEvent,
      nvars = e_mod_stuff$K,
      default_scale = 2.5,
      link = NULL,
      ok_dists = ok_dists
    )  
  
  e_user_prior_intercept_stuff <- e_prior_intercept_stuff <- 
    handle_glm_prior(
      priorEvent_intercept,
      nvars = 1,
      default_scale = 50,
      link = NULL,
      ok_dists = ok_intercept_dists
    )  
  
  e_user_prior_aux_stuff <- e_prior_aux_stuff <-
    handle_glm_prior(
      priorEvent_aux,
      nvars = basehaz$df,
      default_scale = 50,
      link = NULL,
      ok_dists = ok_e_aux_dists
    )
  
  # Priors for association parameters
  a_user_prior_stuff <- a_prior_stuff <- 
    handle_glm_prior(
      priorAssoc,
      nvars = a_K,
      default_scale = 2.5,
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
  if (a_K)
    a_prior_stuff <- autoscale_prior(a_prior_stuff, a_mod_stuff, QR = QR, use_x = FALSE, 
                                     assoc = assoc, family = family)
  
  #-------------------------
  # Data for export to Stan
  #-------------------------
  
  #----- Priors
  
  standata <- list(
    prior_dist                = fetch_array(y_prior_stuff, "prior_dist"), 
    prior_dist_for_intercept  = fetch_array(y_prior_intercept_stuff, "prior_dist"),  
    prior_dist_for_aux        = fetch_array(y_prior_aux_stuff, "prior_dist"),
    e_prior_dist              = e_prior_stuff$prior_dist,
    e_prior_dist_for_intercept= e_prior_intercept_stuff$prior_dist,
    e_prior_dist_for_aux      = e_prior_aux_stuff$prior_dist,
    a_prior_dist              = a_prior_stuff$prior_dist,    
 
    # hyperparameters for longitudinal submodel priors
    prior_mean                 = fetch_array(y_prior_stuff,           "prior_mean"), 
    prior_scale                = fetch_array(y_prior_stuff,           "prior_scale"), 
    prior_df                   = fetch_array(y_prior_stuff,           "prior_df"), 
    prior_mean_for_intercept   = fetch_array(y_prior_intercept_stuff, "prior_mean"),
    prior_scale_for_intercept  = fetch_array(y_prior_intercept_stuff, "prior_scale"), 
    prior_df_for_intercept     = fetch_array(y_prior_intercept_stuff, "prior_df"),  
    prior_mean_for_aux         = fetch_array(y_prior_aux_stuff,       "prior_mean"),
    prior_scale_for_aux        = fetch_array(y_prior_aux_stuff,       "prior_scale"),
    prior_df_for_aux           = fetch_array(y_prior_aux_stuff,       "prior_df"),
    global_prior_scale         = fetch_array(y_prior_stuff, "global_prior_scale"),
    global_prior_df            = fetch_array(y_prior_stuff, "global_prior_df"), 
    
    # hyperparameters for event submodel priors
    e_prior_mean               = e_prior_stuff$prior_mean, 
    e_prior_scale              = e_prior_stuff$prior_scale, 
    e_prior_df                 = e_prior_stuff$prior_df, 
    e_prior_mean_for_intercept = c(e_prior_intercept_stuff$prior_mean),
    e_prior_scale_for_intercept= c(e_prior_intercept_stuff$prior_scale), 
    e_prior_df_for_intercept   = c(e_prior_intercept_stuff$prior_df),
    e_prior_mean_for_aux       = if (basehaz$type == 1L) as.array(0) else 
                                   as.array(e_prior_aux_stuff$prior_mean),  
    e_prior_scale_for_aux      = e_prior_aux_stuff$prior_scale, 
    e_prior_df_for_aux         = e_prior_aux_stuff$prior_df,
    e_global_prior_scale       = e_prior_stuff$global_prior_scale,
    e_global_prior_df          = e_prior_stuff$global_prior_df,
    
    # hyperparameters for assoc parameter priors
    a_prior_mean               = a_prior_stuff$prior_mean, 
    a_prior_scale              = a_prior_stuff$prior_scale, 
    a_prior_df                 = a_prior_stuff$prior_df, 
    a_global_prior_scale       = a_prior_stuff$global_prior_scale,
    a_global_prior_df          = a_prior_stuff$global_prior_df,
    
    # flags
    prior_PD = as.integer(prior_PD)
  )
  
  # prior flag (same prior for all long submodel)
  standata$prior_special_case <- as.integer(
    (length(unique(standata$prior_dist)) == 1L) && all(standata$prior_dist %in% c(0,1,2)))
  
  #----- Longitudinal submodels
  
  # Dimensions and other stuff
  standata$M      <- as.integer(M)
  standata$KM     <- fetch_array(y_mod_stuff, "K")
  standata$NM     <- fetch_array(y_mod_stuff, "N") 
  standata$NM_real<- fetch_array(y_mod_stuff, "real_N") 
  standata$NM_int <- fetch_array(y_mod_stuff, "int_N") 
  standata$K      <- as.integer(sum(standata$KM))
  standata$e_K    <- as.integer(e_mod_stuff$K)
  standata$a_K    <- as.integer(a_K)  
  standata$N      <- as.integer(sum(standata$NM))  
  standata$N_real <- as.integer(sum(fetch_(y_mod_stuff, "real_N"))) 
  standata$N_int  <- as.integer(sum(fetch_(y_mod_stuff, "int_N")))
  standata$N01    <- as.array(t(sapply(fetch(y_mod_stuff, "N01"), cbind))) 
  standata$has_weights       <- as.integer(!is.null(weights))
  standata$has_offset        <- as.integer(!is.null(offset))
  standata$has_intercept     <- fetch_array(y_mod_stuff, "has_intercept")
  standata$has_intercept_nob <- fetch_array(y_mod_stuff, "has_intercept_unbound")
  standata$has_intercept_lob <- fetch_array(y_mod_stuff, "has_intercept_lobound")
  standata$has_intercept_upb <- fetch_array(y_mod_stuff, "has_intercept_upbound")
  standata$has_aux           <- fetch_array(y_mod_stuff, "has_aux")
  standata$xbar              <- fetch_array(y_mod_stuff, "xbar")
  standata$weights <- 
    if (!is.null(weights)) as.array(unlist(y_weights)) else as.array(numeric(0))
  standata$offset  <- 
    if (!is.null(offset)) stop("bug found. offset not yet implemented.") else double(0)
  standata$link    <- as.array(link)
  standata$dense_X <- !sparse
  standata$special_case <- as.integer(FALSE)
  
  # Not used
  standata$K_smooth   <- 0L
  standata$S          <- matrix(NA_real_, standata$N, 0L)
  standata$smooth_map <- integer(0)  
  
  # Design matrices
  X <- as.matrix(Matrix::bdiag(fetch(y_mod_stuff, "xtemp")))
  if (sparse) {
    parts <- extract_sparse_parts(X)
    standata$nnz_X <- length(parts$w)
    standata$w_X <- parts$w
    standata$v_X <- parts$v
    standata$u_X <- parts$u
    standata$X <- array(0, dim = c(0L, dim(X)))
  } else {
    standata$X <- array(X, dim = c(1L, dim(X)))
    standata$nnz_X <- 0L
    standata$w_X <- double(0)
    standata$v_X <- integer(0)
    standata$u_X <- integer(0)
  }  
  
  # Combined response vector
  y_y <- fetch(y_mod_stuff, "y")
  y_is_real <- fetch_(y_mod_stuff, "is_real")
  standata$y_real <- as.array(as.numeric(unlist(y_y[y_is_real])))
  standata$y_int  <- as.array(as.integer(unlist(y_y[!y_is_real])))
  
  # Indexing for combined beta vector, response vector, design matrix, weights, etc
  standata$idx      <- get_idx_array(standata$NM)
  standata$idx_real <- get_idx_array(standata$NM_real)
  standata$idx_int  <- get_idx_array(standata$NM_int)
  standata$idx_K    <- get_idx_array(standata$KM)
  
  # Sum dimensions
  for (i in c("has_aux", paste0("has_intercept", c("", "_nob", "_lob", "_upb")))) 
    standata[[paste0("sum_", i)]] <- as.integer(sum(standata[[i]]))
  
  # Data for group-specific terms
  group <- lapply(y_mod_stuff, function(x) {
    pad_reTrms(Ztlist = x$Ztlist, 
               cnms   = x$cnms, 
               flist  = x$flist)})
  Z              <- fetch(group, "Z")
  y_cnms         <- fetch(group, "cnms")
  y_flist_padded <- fetch(group, "flist")
  t <- length(cnms_nms)   # num of unique grouping factors
  t_i <- which(cnms_nms == id_var) # index of patient-level grouping factor
  pmat <- matrix(0, t, M) # num of group-specific terms
  lmat <- matrix(0, t, M) # num of factor levels
  l <- c()
  for (i in 1:t) {
    for (j in 1:M) {
      pmat[i,j] <- length(y_cnms[[j]][[cnms_nms[i]]])
      lmat[i,j] <- nlevels(y_flist_padded[[j]][[cnms_nms[i]]])
    }
    l[i] <- max(lmat[i,])
    if (!all(lmat[i,] %in% c(0, l[i])))
      stop("The number of factor levels for each of the grouping factors ",
           "must be the same in each of the longitudinal submodels.")     
  }
  qmat <- l * pmat
  p  <- rowSums(pmat) # num group-specific terms for each grouping factor 
  q1 <- rowSums(qmat) # num group-specific coefs for each grouping factor
  q2 <- colSums(qmat) # num group-specific coefs for each submodel
  q  <- sum(qmat)     # total num group-specific coefs
  b_nms <- unlist(Map(make_b_nms, group, m = seq(M)))
  g_nms <- unlist(
    lapply(1:M, FUN = function(m) {
      lapply(1:length(group[[m]]$cnms), FUN = function(i) {
        paste(paste0("Long", m), group[[m]]$cnms[[i]], names(group[[m]]$cnms)[i], sep = "|")
      })
    })
  )
  standata$t    <- as.integer(t)
  standata$t_i  <- as.integer(t_i)
  standata$pmat <- as.array(pmat)
  standata$p    <- as.array(p)
  standata$l    <- as.array(l)
  standata$qmat <- as.array(qmat)
  standata$q1   <- as.array(q1)
  standata$q2   <- as.array(q2)
  standata$q    <- as.integer(q)
  standata$len_theta_L <- sum(choose(p, 2), p)
  Zmerge <- Matrix::bdiag(Z)
  parts <- rstan::extract_sparse_parts(Zmerge)
  standata$num_non_zero <- as.integer(length(parts$w))
  standata$w <- parts$w
  standata$v <- parts$v
  standata$u <- as.array(parts$u)
  
  # Hyperparameters for decov prior
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
  
  # Families
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
  
  #----- Event submodel (including GK quadrature)
  
  # Dimensions, response, design matrix, etc
  standata$Npat            <- as.integer(e_mod_stuff$Npat)
  standata$quadnodes       <- as.integer(quadnodes)
  standata$quadweight      <- as.array(e_mod_stuff$quadweight)
  standata$Npat_times_quadnodes <- as.integer(e_mod_stuff$Npat * quadnodes)
  standata$nrow_y_Xq       <- NROW(a_mod_stuff[[1]]$mod_eta$xtemp)
  standata$nrow_e_Xq       <- NROW(e_mod_stuff$xtemp)
  standata$e_times         <- c(e_mod_stuff$eventtime, unlist(e_mod_stuff$quadpoint))
  standata$e_d             <- c(e_mod_stuff$d, rep(1, length(unlist(e_mod_stuff$quadpoint))))
  standata$e_has_intercept <- as.integer(e_has_intercept)
  standata$e_Xq            <- e_mod_stuff$xtemp
  standata$e_xbar          <- as.array(e_mod_stuff$xbar)
  standata$e_weights       <- as.array(e_weights)
  standata$e_weights_rep   <- as.array(rep(e_weights, times = quadnodes))
  
  # Baseline hazard
  standata$basehaz_type <- as.integer(basehaz$type)
  standata$basehaz_df   <- as.integer(basehaz$df)
  if (basehaz$type_name == "weibull") {
    standata$basehaz_X <- matrix(log(standata$e_times), length(standata$e_times), 1) 
  } else if (basehaz$type_name == "bs") {
    standata$basehaz_X <- as.array(predict(basehaz$bs_basis, standata$e_times)) 
  } else if (basehaz$type_name == "piecewise") {
    e_times_quantiles <- cut(standata$e_times, basehaz$knots, 
                             include.lowest = TRUE, labels = FALSE)
    tmp <- matrix(NA, length(e_times_quantiles), basehaz$df)
    for (i in 1:basehaz$df) 
      tmp[, i] <- ifelse(e_times_quantiles == i, 1, 0)
    standata$basehaz_X <- as.array(tmp)
  } else {
    standata$basehaz_X <- matrix(0,0,0)  
  }
  
  #----- Association structure
  
  standata$assoc <- as.integer(a_K > 0L) # any association structure, 1 = yes

  # Indicator for which components are required to build the association terms
  standata$assoc_uses <- sapply(
    c("etavalue", "etaslope", "etaauc", "muvalue", "muslope", "muauc"), 
    function(x, assoc) {
      nm_check <- switch(x,
                         etavalue = "^eta|^mu",
                         etaslope = "etaslope|muslope",
                         etaauc   = "etaauc|muauc",
                         muvalue  = "muvalue|muslope",
                         muslope  = "muslope",
                         muauc    = "muauc")
      sel <- grep(nm_check, rownames(assoc))
      as.integer(any(unlist(assoc[sel,])))
    }, assoc = assoc)  
    
  # Indexing for desired association types
  # !!! must be careful with corresponding use of indexing in Stan code
  # 1 = ev; 2 = es; 3 = ea; 4 = mv; 5 = ms; 6 = ma;
  # 7 = shared_b; 8 = shared_coef;
  # 9 = ev_data; 10 = es_data; 11 = mv_data; 12 = ms_data;
  # 13 = evev; 14 = evmv; 15 = mvev; 16 = mvmv;
  sel <- grep("which|null", rownames(assoc), invert = TRUE)
  standata$has_assoc <- matrix(as.integer(assoc[sel,]), ncol = M) 
  
  # Data for calculating value, slope, auc in GK quadrature 
  for (i in c("eta", "eps", "auc")) {
    nm_check <- switch(i,
                       eta = "^eta|^mu",
                       eps = "slope",
                       auc = "auc")
    sel <- grep(nm_check, rownames(assoc))
    if (any(unlist(assoc[sel,]))) {
      tmp_stuff <- fetch(a_mod_stuff, paste0("mod_", i))
      X_tmp <- as.matrix(Matrix::bdiag(fetch(tmp_stuff, "xtemp")))
      group_tmp <- lapply(tmp_stuff, function(x) {
        pad_reTrms(Ztlist = x[["group"]][["Ztlist"]], 
                   cnms   = x[["group"]][["cnms"]], 
                   flist  = x[["group"]][["flist"]])})
      Z_tmp <- Matrix::bdiag(fetch(group_tmp, "Z"))      
    } else {
      X_tmp <- matrix(0,0,standata$K)
      Z_tmp <- matrix(0,0,0) 
    }
    parts_Z_tmp <- rstan::extract_sparse_parts(Z_tmp)
    standata[[paste0("y_Xq_", i)]] <- as.array(X_tmp)
    standata[[paste0("nnz_Zq_", i)]] <- as.integer(length(parts_Z_tmp$w))
    standata[[paste0("w_Zq_", i)]] <- parts_Z_tmp$w
    standata[[paste0("v_Zq_", i)]] <- parts_Z_tmp$v
    standata[[paste0("u_Zq_", i)]] <- as.array(parts_Z_tmp$u)    
  }  
  
  # Data for slope association structure
  standata$eps <- eps
  
  # Data for auc association structure
  standata$auc_quadnodes <- 
    as.integer(auc_quadnodes)
  standata$Npat_times_auc_quadnodes <- 
    as.integer(e_mod_stuff$Npat * auc_quadnodes)  
  standata$nrow_y_Xq_auc <- 
    as.integer(NROW(a_mod_stuff[[1]]$mod_auc$xtemp))
  auc_quadweights <- 
    unlist(lapply(e_mod_stuff$quadtimes, function(x) 
      lapply(x, function(y) 
        lapply(get_quadpoints(auc_quadnodes)$weights, unstandardise_quadweights, 0, y)))) 
  standata$auc_quadweights <- 
    if (standata$assoc_uses[3]) as.array(auc_quadweights) else double(0)

  # Interactions between association terms and data
  # design matrix for the interactions
  standata$y_Xq_data <- 
    as.array(t(as.matrix(do.call("cbind", fetch(a_mod_stuff, "xmat_data")))))
  # number of columns in y_Xq_data corresponding to each interaction type 
  # (ie, etavalue, etaslope, muvalue, muslope) for each submodel
  standata$a_K_data  <- fetch_array(a_mod_stuff, "K_data")  
  
  # Interactions between association terms
  standata$which_interactions      <- as.array(unlist(assoc["which_interactions",]))
  standata$size_which_interactions <- c(sapply(assoc["which_interactions",], sapply, length))
   
  # Shared random effects
  standata$which_b_zindex    <- as.array(unlist(assoc["which_b_zindex",]))
  standata$which_coef_zindex <- as.array(unlist(assoc["which_coef_zindex",]))
  standata$which_coef_xindex <- as.array(unlist(assoc["which_coef_xindex",]))
  standata$size_which_b      <- as.array(sapply(assoc["which_b_zindex",    ], length))
  standata$size_which_coef   <- as.array(sapply(assoc["which_coef_zindex", ], length))

  # Sum dimensions
  for (i in c("a_K_data", paste0("size_which_", c("b", "coef", "interactions")))) {
    standata[[paste0("sum_", i)]] <- as.integer(sum(standata[[i]]))
  }
  
  #----------------
  # Initial values
  #----------------
  
  if (is.character(init) && (init =="prefit")) {
    init <- generate_init_function(y_mod_stuff, e_mod_stuff, standata)
  } else if (is.character(init) && (init == "prefit_vb")) {
    cat("Obtaining initial values using variational bayes\n")
    dropargs <- c("chains", "cores", "iter")
    vbdots <- list(...)
    for (i in dropargs) 
      vbdots[[i]] <- NULL
    vbargs <- c(list(stanmodels$mvmer, pars = "mean_PPD", data = standata, 
                     algorithm = "meanfield"), vbdots)
    initfit <- do.call(rstan::vb, vbargs)
    initmeans <- rstan::get_posterior_mean(initfit)
    initnms <- rownames(initmeans)
    inits <- generate_init_function(y_mod_stuff, e_mod_stuff, standata)()
    sel <- c("gamma_nob", "gamma_lob", "gamma_upb", "z_beta", "aux_unscaled", 
             "z_b", "z_T", "rho", "zeta", "tau", "global", "local2", "local4", 
             "mix", "ool", "noise")
    for (i in sel) {
      sel_i <- grep(paste0("^", i, "\\."), initnms)
      if (length(sel_i))
        inits[[i]] <- as.array(initmeans[sel_i,])
    }
    init <- function() inits
  }
  
  #---------------
  # Prior summary
  #---------------
  
  prior_info <- summarize_jm_prior(
    user_priorLong = y_user_prior_stuff,
    user_priorLong_intercept = y_user_prior_intercept_stuff,
    user_priorLong_aux = y_user_prior_aux_stuff,
    user_priorEvent = e_user_prior_stuff,
    user_priorEvent_intercept = e_user_prior_intercept_stuff,
    user_priorEvent_aux = e_user_prior_aux_stuff,
    user_priorAssoc = a_user_prior_stuff,
    user_prior_covariance = prior_covariance,
    y_has_intercept = sapply(y_mod_stuff, `[[`, "has_intercept"),
    y_has_predictors = sapply(y_mod_stuff, function(x) x$K > 0),
    e_has_intercept = e_mod_stuff$has_intercept,
    e_has_predictors = e_mod_stuff$K > 0,
    has_assoc = a_K > 0,
    adjusted_priorLong_scale = fetch(y_prior_stuff, "prior_scale"),
    adjusted_priorLong_intercept_scale = fetch(y_prior_intercept_stuff, "prior_scale"),
    adjusted_priorLong_aux_scale = fetch(y_prior_aux_stuff, "prior_scale"),
    adjusted_priorEvent_scale = e_prior_stuff$prior_scale,
    adjusted_priorEvent_intercept_scale = e_prior_intercept_stuff$prior_scale,
    adjusted_priorEvent_aux_scale = e_prior_aux_stuff$prior_scale,
    adjusted_priorAssoc_scale = a_prior_stuff$prior_scale,
    family = family, basehaz = basehaz
  )  

  #-----------
  # Fit model
  #-----------
  
  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$jm
  pars <- c(if (standata$sum_has_intercept) "alpha",
            if (standata$K) "beta",
            if (standata$e_has_intercept) "e_gamma",
            if (standata$e_K) "e_beta",
            if (standata$a_K) "a_beta",
            if (standata$q) "b",
            if (standata$sum_has_aux) "aux",
            if (length(standata$basehaz_X)) "e_aux",
            if (standata$len_theta_L) "theta_L",
            "mean_PPD")
            
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
  check_stanfit(stanfit)
  
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
    if (assoc["etaauc",           ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|etaauc"))
    if (assoc["muvalue",          ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue"))
    if (assoc["muvalue_data",     ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue:", colnames(a_mod_stuff[[m]][["xq_data"]][["muvalue_data"]])))    
    if (assoc["muvalue_etavalue", ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_etavalue"]], "|etavalue"))
    if (assoc["muvalue_muvalue",  ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muvalue:Long", assoc["which_interactions",][[m]][["muvalue_muvalue"]],  "|muvalue"))
    if (assoc["muslope",          ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muslope"))
    if (assoc["muslope_data",     ][[m]]) a_nms <- c(a_nms, paste0("Assoc|Long", m,"|muslope:", colnames(a_mod_stuff[[m]][["xq_data"]][["muslope_data"]])))    
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
  
  # Sigma values in stanmat, and Sigma names
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
                 if (length(standata$q)) c(paste0("b[", b_nms, "]")),
                 y_aux_nms,
                 e_aux_nms,
                 if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                 paste0("Long", 1:M, "|mean_PPD"), 
                 "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  
  n_grps <- standata$l - 1
  names(n_grps) <- cnms_nms  # n_grps is num. of levels within each grouping factor
  names(p) <- cnms_nms       # p is num. of variables within each grouping factor
  
  # Undo ordering of matrices if bernoulli
  y_mod_stuff <- lapply(y_mod_stuff, unorder_bernoulli)
  
  fit <- nlist(stanfit, family, formula = list_nms(c(formulaLong, formulaEvent), M), 
               id_var, time_var, offset, weights, quadnodes, basehaz,
               M, cnms, n_yobs = unlist(list_nms(fetch(y_mod_stuff, "N"), M)), 
               n_subjects = e_mod_stuff$Npat, n_grps, assoc,
               y_mod_stuff, e_mod_stuff, a_mod_stuff,
               fr = list_nms(c(fetch(a_mod_stuff, "model_frame"), 
                               list(e_mod_stuff$model_frame)), M),
               y = list_nms(fetch(y_mod_stuff, "y"), M),
               d = e_mod_stuff$d, eventtime = e_mod_stuff$eventtime,
               epsilon = if (standata$assoc_uses[2]) eps else NULL,
               dataLong, dataEvent, call, na.action, algorithm, 
               standata = NULL, terms = NULL, prior.info = prior_info,
               stan_function = "stan_jm")
  out <- stanmvreg(fit)
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
handle_glmod <- function(mc, family, supported_families, supported_links, 
                         sparse, env = parent.frame()) {
  
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
  mod <- suppressWarnings(eval(mc, envir = env))     	
  
  # Response vector
  y <- as.vector(lme4::getME(mod, "y"))
  y <- validate_glm_outcome_support(y, family)
  if (is.binomial(family$family) && NCOL(y) == 2L)
    STOP_binomial()
  
  # Design matrix
  xtemp <- xbar <- has_intercept <- NULL # useless assignments to pass R CMD check
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
  group <- eval(group_call, envir = env)$reTrms      	
  
  Ztlist <- group$Ztlist
  cnms   <- group$cnms
  flist  <- group$flist
  
  # Reorder y, X, Z if bernoulli (zeros first)
  if (is.binomial(family$family) && all(y %in% 0:1)) {      
    ord    <- order(y)
    y      <- y     [ord]
    xtemp  <- xtemp [ord, , drop = FALSE]  
    Ztlist <- lapply(Ztlist, function(x) x[, ord, drop = FALSE]) 
    flist  <- flist[ord, , drop = FALSE] 
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
  if (is.binomial(famname) && !is_bernoulli)
    STOP_binomial()
  
  # Require intercept for certain family and link combinations
  if (!has_intercept) {
    link <- which(supported_links == family$link)
    linkname <- supported_links[link]
    needs_intercept <- !is_gaussian && linkname == "identity" ||
      is_gamma && linkname == "inverse" ||
      is.binomial(famname) && linkname == "log"
    if (needs_intercept)
      stop("To use the specified combination of family and link (", famname, 
           ", ", linkname, ") the model must have an intercept.")
  }    
    
  # Return list
  nlist(mod, is_real, y, x, xtemp, xbar, N01, ord, famname,
    offset, Ztlist, cnms, flist, has_intercept, has_intercept_unbound,
    has_intercept_lobound, has_intercept_upbound, has_aux, N, real_N, int_N, K,
    is_bernoulli, is_nb, is_gaussian, is_gamma, is_ig, is_continuous, is_lmer)
}

# Function to return a single cnms object for all longitudinal submodels
#
# @param x A list, with each element being a cnms object returned by (g)lmer
get_common_cnms <- function(x, stub = "Long") {
  nms <- lapply(x, names)
  unique_nms <- unique(unlist(nms))
  cnms <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    unlist(lapply(1:length(x), function(m) 
      if (nm %in% nms[[m]]) paste0(stub, m, "|", x[[m]][[nm]])))
  })
  names(cnms) <- unique_nms
  cnms
}

# Function to substitute variables in the formula of a fitted model
# with the corresponding predvars based on the terms object for the model.
# (This is useful since lme4::glFormula doesn't allow a terms object to be 
# passed as the first argument instead of a model formula).
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
    fm <- reformulate(fm[[3L]], response = formula(mod)[[2L]])
  }
  fm
}

# Check the id_var argument is valid and is included appropriately in the
# formulas for each of the longitudinal submodels
#
# @param id_var The character string that the user specified for the id_var
#   argument -- will have been set to NULL if the argument was missing.
# @param y_cnms A list of length M with the cnms for each longitudinal submodel
# @param y_flist A list of length M with the flist for each longitudinal submodel
# @return Returns the character string corresponding to the appropriate id_var.
#   This will either be the user specified id_var argument or the only grouping
#   factor.
check_id_var <- function(id_var, y_cnms, y_flist) {
  len_cnms <- sapply(y_cnms, length)
  if (any(len_cnms > 1L)) {  # more than one grouping factor
    if (is.null(id_var)) {
      stop("'id_var' must be specified when using more than one grouping factor",
           call. = FALSE)
    } else {
      lapply(y_cnms, function(x)  if (!(id_var %in% names(x)))
        stop("'id_var' must be included as a grouping factor in each ",
             "of the longitudinal submodels", call. = FALSE)) 
      mapply(function(cnms, flist, id_var) { # loop over submodels
        # If there is more than one grouping factor, then make sure that the id_var 
        # is the lowest level of clustering, by checking that the observations for 
        # a given individual don't have more than one level for each other grouping
        # factor included in the model
        nms <- grep(id_var, names(cnms), value = TRUE, invert = TRUE)
        if (length(nms)) { # submodel has additional grouping factors
          lapply(nms, function(x) { # loop over additional grouping factors
            # within each ID, count the number of levels for the additional grouping factor 
            tally <- tapply(flist[[x]], flist[[id_var]], function(y) length(unique(y)))
            # within each ID, ensure max of 1 level for each other grouping factor
            if (!all(tally == 1L))
              stop("The 'id_var' must correspond to the lowest level of clustering. ",
                   "Yet for some levels of '", id_var, "' there appears to be more ",
                   "than one level for '", x, "'.", call. = FALSE)
          }) 
        }
      }, cnms = y_cnms, flist = y_flist, MoreArgs = list(id_var = id_var))    
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
    mod_stuff$weights <- mod_stuff$weights[order(mod_stuff$ord)]
    mod_stuff$xtemp   <- mod_stuff$xtemp[order(mod_stuff$ord), , drop = FALSE]
    mod_stuff$Ztlist  <- lapply(mod_stuff$Ztlist, function(x) x[, order(mod_stuff$ord), drop = FALSE])
    mod_stuff$flist   <- mod_stuff$flist[order(mod_stuff$ord), , drop = FALSE]
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
handle_coxmod <- function(mc, quadnodes, id_var, unique_id_list, sparse,
                          env = parent.frame()) {
  
  mc[[1]] <- quote(survival::coxph) 
  mc$x    <- TRUE
  mod <- eval(mc, envir = env)
  mf1 <- expand.model.frame(mod, id_var, na.expand = TRUE)
  # since lme4 promotes character grouping variables to factors
  if (is.character(mf1[[id_var]]))
    mf1[[id_var]] <- as.factor(mf1[[id_var]])
  mf2 <- cbind(unclass(mf1[,1]), mf1[, -1, drop = FALSE])
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
  xtemp <- xbar <- has_intercept <- NULL # useless assignments for R CMD check
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
                           ok_assoc_interactions, lag, id_var, M) {

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
      STOP_combination_not_allowed(assoc, "etaauc",   "muauc")
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
  
  # Lag for association structure
  assoc$which_lag <- lag

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
  x <- gsub("^shared_b\\(.*",    "shared_b",    x) 
  x <- gsub("^shared_coef\\(.*", "shared_coef", x) 
  for (i in ok_assoc_data)
    x <- gsub(paste0("^", i, "_data\\(.*"),    paste0(i, "_data"), x)
  for (i in ok_assoc_interactions) for (j in ok_assoc_interactions)
    x <- gsub(paste0("^", i, "_", j, "\\(.*"), paste0(i, "_", j),  x) 
  x     
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
#   on the desired association structure for all longitudinal submodels)
# @param id_var The name on the ID variable
# @param time_var The name of the time variable
# @param eps The time shift used for the numerical derivative calculation 
#   based on a one-sided different
# @param dataAssoc An optional data frame containing data for interactions within
#   the association terms
handle_assocmod <- function(m, mc, dataLong, y_mod_stuff, id_list, times, assoc, 
                            id_var, time_var, eps, auc_quadnodes, 
                            dataAssoc = NULL, env = parent.frame()) {
  
  # Obtain a model frame defined as a data.table
  rows <- rownames(model.frame(y_mod_stuff$mod))
  df   <- as.data.frame(dataLong)[rows,]
  mf   <- data.table::data.table(df, key = c(id_var, time_var))
  mf[[time_var]] <- as.numeric(mf[[time_var]]) # ensure no rounding on merge
  
  # Update longitudinal submodel formula to reflect predvars
  mc$formula <- use_predvars(y_mod_stuff$mod)

  # Design matrices for calculating eta, eps, lag, auc and data interactions
  # in association structure
  mc[[1]] <- quote(lme4::glFormula)
  parts <- make_assoc_parts(newdata = mf, assoc = assoc, m = m, id_var = id_var, 
                            time_var = time_var, id_list = id_list, times = times, 
                            eps = eps, auc_quadnodes = auc_quadnodes,
                            use_function = handle_glFormula, 
                            mc = mc, y_mod_stuff = y_mod_stuff, 
                            dataAssoc = dataAssoc, env = env)
  
  # If association structure is based on shared random effects or shared 
  # coefficients then construct a matrix with the estimated b parameters
  # from the separate glmod (for the id_var grouping factor only). Note this
  # matrix is not passed to standata, but just used for autoscaling the 
  # priors for association parameters.
  sel_shared <- grep("^shared", rownames(assoc))
  if (any(unlist(assoc[sel_shared,]))) {
    # flist for long submodel
    flist_tmp <- lme4::getME(y_mod_stuff$mod, "flist")
    # which grouping factor is id_var
    Gp_sel <- which(names(flist_tmp) == id_var) 
    # grouping factor indices
    Gp <- lme4::getME(y_mod_stuff$mod, "Gp")  
    b_beg <- Gp[[Gp_sel]] + 1
    b_end <- Gp[[Gp_sel + 1]]
    # b vector for grouping factor = id_var
    b_vec <- lme4::getME(y_mod_stuff$mod, "b")[b_beg:b_end]
    # convert to Npat * n_re matrix
    b_mat <- matrix(b_vec, nrow = length(levels(flist_tmp[[Gp_sel]])), byrow = TRUE)
  } else b_mat <- NULL
  
  parts$b_mat <- b_mat
  return(parts)
}

# Function to construct quantities, primarily design matrices (X, Zt), that
# will be used to evaluate the longitudinal submodel contributions to the 
# association structure in the event submodel. For example, the design matrices
# evaluated at the quadpoints, the quadpoints, lagged quadpoints, auc quadpoints,
# and so on. Exactly what quantities are returned depends on what is specified
# in the use_function argument.
#
# @param data A model frame used for constructing the design matrices
# @param assoc A named list returned by a call to validate_assoc (details
#   on the desired association structure for all longitudinal submodels)
# @param id_var The name on the ID variable
# @param time_var The name of the time variable
# @param id_list A vector of subject IDs
# @param times A vector (or possibly a list of vectors) of times at which the 
#   design matrices should be evaluated (most likely the event times and the
#   quadrature times)
# @param eps A numeric value used as the time shift for numerically evaluating
#   the slope of the longitudinal submodel using a one-sided difference
# @param auc_quadnodes An integer specifying the number of quadrature nodes to
#   use when evaluating the area under the curve for the longitudinal submodel
# @param use_function The function to call which will return the design 
#   matrices for eta, eps, lag, auc, etc.
# @param ... Additional arguments passes to use_function
# @return A named list
make_assoc_parts <- function(newdata, assoc, id_var, time_var, 
                             id_list, times, eps = 1E-5, auc_quadnodes = 15L, 
                             dataAssoc = NULL, use_function = handle_glFormula, 
                             ...) {
  dots <- list(...)
  m <- dots$m 
  if (is.null(m)) stop("Argument m must be specified in dots.")

  # Apply lag
  lag <- assoc["which_lag",][[m]]
  if (!lag == 0) {
    times <- lapply(times, function(x, lag) {
      newtimes <- x - lag
      newtimes[newtimes < 0] <- 0.0  # use baseline where lagged t is before baseline
      newtimes
    }, lag = lag)  
  }
  
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
  dataQ <- rolling_merge(data = newdata, ids = id_list, times = times)
  mod_eta <- use_function(newdata = dataQ, ...)
  
  # If association structure is based on slope, then calculate design 
  # matrices under a time shift of epsilon
  sel_slope <- grep("etaslope|muslope", rownames(assoc))
  if (any(unlist(assoc[sel_slope,]))) {
    dataQ_eps <- dataQ
    dataQ_eps[[time_var]] <- dataQ_eps[[time_var]] + eps
    mod_eps <- use_function(newdata = dataQ_eps, ...)
  } else mod_eps <- NULL 
  
  # If association structure is based on area under the marker trajectory, then 
  # calculate design matrices at the subquadrature points
  sel_auc <- grep("etaauc|muauc", rownames(assoc))
  if (any(unlist(assoc[sel_auc,]))) {
    # Return a design matrix that is (quadnodes * auc_quadnodes * Npat) rows 
    auc_quadtimes <- 
      lapply(times, function(x) unlist(
        lapply(x, function(y) 
          lapply(get_quadpoints(auc_quadnodes)$points, unstandardise_quadpoints, 0, y))))
    auc_quadweights <- 
      lapply(times, function(x) unlist(
        lapply(x, function(y) 
          lapply(get_quadpoints(auc_quadnodes)$weights, unstandardise_quadweights, 0, y))))
    ids2 <- rep(id_list, each = auc_quadnodes)
    dataQ_auc <- rolling_merge(data = newdata, ids = ids2, times = auc_quadtimes)
    mod_auc <- use_function(newdata = dataQ_auc, ...)
  } else mod_auc <- auc_quadtimes <- auc_quadweights <- NULL
  
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
                        if ("y_mod_stuff" %in% names(dots)) { # call from stan_jm
                          y_mod_stuff <- dots$y_mod_stuff
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
                        } else { # call from posterior_survfit
                          object <- dots$object
                          oldcall <- getCall(object$glmod_stuff[[m]]$mod)
                          naa <- oldcall$na.action
                          subset <- eval(oldcall$subset) 
                          if (is.null(dataAssoc)) {
                            df <- 
                              tryCatch(eval(oldcall$data), error = function(e) {
                                tryCatch(if (is(object$dataLong, "list")) 
                                  object$dataLong[[m]] else object$dataLong, error = function(e) {
                                    stop("Bug found: cannot find data frame to use for ",
                                         "assoc data interactions, please report bug.")
                                  })
                              })
                          } else df <- dataAssoc
                          mf2 <- eval(call("model.frame", ff, data = df, subset = subset, 
                                           na.action = naa))                          
                        }  
                        mf2 <- data.table::data.table(mf2, key = c(id_var, time_var))
                        mf2[[time_var]] <- as.numeric(mf2[[time_var]])
                        mf2q <- rolling_merge(data = mf2, ids = id_list, times = times)
                        xq <- stats::model.matrix(fm, data = mf2q)
                        if ("(Intercept)" %in% colnames(xq)) xq <- xq[, -1L, drop = FALSE]
                        if (!ncol(xq))
                          stop(paste0("Bug found: A formula was specified for the '", x, "' association ", 
                                      "structure, but the resulting design matrix has no columns."), call. = FALSE)
                      } else xq <- matrix(0, length(unlist(times)), 0)
                      xq
                    }, simplify = FALSE, USE.NAMES = TRUE)
  K_data <- sapply(xq_data, ncol)
  xmat_data <- do.call(cbind, xq_data)
  
  ret <- nlist(times, mod_eta, mod_eps, mod_auc, xq_data, xmat_data, K_data)
  
  structure(ret, times = times, lag = lag, eps = eps, auc_quadnodes = auc_quadnodes,
            auc_quadtimes = auc_quadtimes, auc_quadweights = auc_quadweights)
}                              

# Carry out a rolling merge
#
# @param data A data.table with a set key corresponding to ids and times
# @param ids A vector of ids to merge against
# @param times A vector of (new) times to merge against
# @return A data.table formed by a merge of ids, times, and the closest 
#   preceding (in terms of times) rows in data
rolling_merge <- function(data, ids, times) {
  if (is(times, "list")) {
    return(do.call(rbind, lapply(times, FUN = function(x) 
      data[list(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)])))      
  } else 
    return(data[list(ids, times), roll = TRUE, rollends = c(TRUE, TRUE)])     
}

# Evaluate a glFormula call and return model components
# 
# @param mc A glFormula call
# @param newdata A data frame to substitute into the data argument of the call
# @param y_mod_stuff A named list, returned by a call to handle_glmod (and
#   containing an indicator of whether the original longitudinal submodel had
#   an intercept term)
# @param m Argument ignored (included to avoid error when passing m to 
#   make_assoc_Xparts function)
# @param env The environment in which to evaluate the glFormula call (note that
#   although the formula and data have been substituted, there may be additional
#   arugments (for example family) that need to be evaluated in the environment
#   of the original stan_jm call).
handle_glFormula <- function(mc, newdata, y_mod_stuff, m = NULL, 
                             env = parent.frame()) { 
  mc$data <- newdata
  mod    <- eval(mc, env)
  x      <- as.matrix(mod$X)
  xtemp  <- if (y_mod_stuff$has_intercept) x[, -1L, drop = FALSE] else x  
  xtemp  <- sweep(xtemp, 2, y_mod_stuff$xbar, FUN = "-")
  group  <- mod$reTrms    
  beta   <- fixef(y_mod_stuff$mod)
  b      <- lme4::getME(y_mod_stuff$mod, "b")
  linpred <- linear_predictor.default(beta, x) # offset not accomodated here
  linpred <- linpred + (t(as.matrix(group$Zt)) %*% b)
  nlist(xtemp, group, linpred)
}   

# Function to calculate the number of association parameters in the model
#
# @param assoc A list of length M with information about the association structure
#   type for each submodel, returned by an mapply call to validate_assoc
# @param a_mod_stuff A list of length M with the design matrices related to
#   the longitudinal submodels in the GK quadrature, returned by an mapply 
#   call to handle_assocmod
# @return Integer indicating the number of association parameters in the model 
get_num_assoc_pars <- function(assoc, a_mod_stuff) {
  sel1 <- c("etavalue", "etaslope", "etaauc", 
            "muvalue", "muslope", "muauc")
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
# @param mod_stuff A named list returned by a call to either handle_glmod,
#   handle_coxmod, or handle_assocmod
# @param QR A logical specifying whether QR decomposition is used for the x matrix
# @param use_x A logical specifying whether to autoscale the priors based on
#   the standard deviations of the predictor variables
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). Cannot
#   be NULL if autoscaling priors for the association parameters.
# @param min_prior_scale The minimum allowed for prior scales
# @return A named list of the same structure as returned by handle_glm_prior
autoscale_prior <- function(prior_stuff, mod_stuff, QR, use_x = FALSE, 
                            min_prior_scale = 1e-12, assoc = NULL, family = NULL) {
  
  is_glmod    <- ("y"         %in% names(mod_stuff))
  is_coxmod   <- ("eventtime" %in% names(mod_stuff))
  is_assocmod <- (!any(is_glmod, is_coxmod))
  
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
  
  if (is_assocmod && prior_stuff$prior_dist > 0L && prior_stuff$prior_autoscale) {
    # Evaluate mean and SD of each of the association terms that will go into
    # the linear predictor for the event submodel (as implicit "covariates").
    # (NB the approximate association terms are calculated using coefs
    # from the separate longitudinal submodels estimated using glmer).
    # The mean will be used for centering each association term.
    # The SD will be used for autoscaling the prior for each association parameter.
    if (is.null(assoc) || is.null(family))
      stop("'assoc' and 'family' cannot be NULL when autoscaling association parameters.")
    assoc_terms <- make_assoc_terms(parts = mod_stuff, assoc = assoc, family = family)
    a_beta_scale <- apply(assoc_terms, 2L, scale_val)
    prior_stuff$prior_scale <- 
      pmax(min_prior_scale, prior_stuff$prior_scale / a_beta_scale) 
  }
  
  prior_stuff$prior_scale <- 
    as.array(pmin(.Machine$double.xmax, prior_stuff$prior_scale))
  
  prior_stuff
}


# Function to construct a design matrix for the association structure in
# the event submodel, to be multiplied by a vector of association parameters
#
# @param assoc An array with information about the desired association 
#   structure, returned by a call to validate_assoc
# @param parts A list equal in length to the number of markers. Each element
#   parts[[m]] should contain a named list with components $mod_eta, $mod_eps,
#   $mod_auc, which each contain either the linear predictor at quadtimes, 
#   quadtimes + eps, and auc quadtimes, or the design matrices
#   used for constructing the linear predictor. Each element parts[[m]] should 
#   also contain $xmat_data and $K_data.
# @param family A list of family objects, equal in length to the number of 
#   longitudinal submodels
# @param ... If parts does not contain the linear predictors, then this should
#   include elements beta and b, each being a length M list of parameters for the
#   longitudinal submodels
# @return A design matrix containing the association terms to be multiplied by
#   the association paramters.
make_assoc_terms <- function(parts, assoc, family, ...) {
  M <- length(parts)
  a_X <- list()
  mark <- 1
  for (m in 1:M) {
    times  <- attr(parts[[m]], "times")
    eps    <- attr(parts[[m]], "eps")  
    qnodes <- attr(parts[[m]], "auc_quadnodes")
    qwts   <- unlist(attr(parts[[m]], "auc_quadweights"))
    
    if (!assoc["null",][[m]]) {
      invlink_m <- family[[m]]$linkinv    
      eta_m <- get_element(parts, m = m, "eta", ...)
      eps_m <- get_element(parts, m = m, "eps", ...)
      auc_m <- get_element(parts, m = m, "auc", ...)
      data_m <- get_element(parts, m = m, "xmat_data", ...)
      K_data_m <- get_element(parts, m = m, "K_data", ...)
      
      # etavalue and any interactions      
      if (assoc["etavalue",][[m]]) { # etavalue
        a_X[[mark]] <- eta_m
        mark <- mark + 1
      }
      if (assoc["etavalue_data",][[m]]) { # etavalue*data
        idx_ev <- which(names(K_data_m) == "etavalue_data")
        cbeg  <- sum(K_data_m[0:(idx_ev-1)]) + 1
        cend  <- sum(K_data_m[0: idx_ev   ])
        val <- as.vector(eta_m) * data_m[, cbeg:cend, drop = FALSE]
        a_X[[mark]] <- val
        mark <- mark + 1
      }
      if (assoc["etavalue_etavalue",][[m]]) { # etavalue*etavalue
        sel <- assoc["which_interactions",][[m]][["etavalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val   <- eta_m * eta_j 
          a_X[[mark]] <- val
          mark <- mark + 1
        }
      } 
      if (assoc["etavalue_muvalue",][[m]]) { # etavalue*muvalue
        sel <- assoc["which_interactions",][[m]][["etavalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- eta_m * invlink_j(eta_j) 
          a_X[[mark]] <- val
          mark <- mark + 1             
        }
      }       
      # etaslope and any interactions
      if (assoc["etaslope",][[m]]) { # etaslope
        dydt_m <- (eps_m - eta_m) / eps
        a_X[[mark]] <- dydt_m
        mark <- mark + 1             
      }
      if (assoc["etaslope_data",][[m]]) { # etaslope*data
        dydt_m <- (eps_m - eta_m) / eps
        idx_es <- which(names(K_data_m) == "etaslope_data")
        cbeg  <- sum(K_data_m[0:(idx_es-1)]) + 1
        cend  <- sum(K_data_m[0: idx_es   ])
        val <- as.vector(dydt_m) * data_m[, cbeg:cend, drop = FALSE]
        a_X[[mark]] <- val
        mark <- mark + 1            
      }
      # etaauc
      if (assoc["etaauc",][[m]]) { # etaauc
        val   <- c()
        for (j in 1:length(eta_m)) {
          wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
          auc_j <- auc_m[((j-1) * qnodes + 1):(j * qnodes)]
          val[j] <- sum(wgt_j * auc_j)
        }
        a_X[[mark]] <- val
        mark <- mark + 1            
      }        
      # muvalue and any interactions
      if (assoc["muvalue",][[m]]) { # muvalue
        val <- invlink_m(eta_m) 
        a_X[[mark]] <- val
        mark <- mark + 1            
      }
      if (assoc["muvalue_data",][[m]]) { # muvalue*data
        idx_mv <- which(names(K_data_m) == "muvalue_data")
        cbeg  <- sum(K_data_m[0:(idx_mv-1)]) + 1
        cend  <- sum(K_data_m[0: idx_mv   ])
        val   <- as.vector(invlink_m(eta_m)) * data_m[, cbeg:cend, drop = FALSE] 
        a_X[[mark]] <- val
        mark <- mark + 1           
      }
      if (assoc["muvalue_etavalue",][[m]]) { # muvalue*etavalue
        sel <- assoc["which_interactions",][[m]][["muvalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val   <- invlink_m(eta_m) * eta_j 
          a_X[[mark]] <- val
          mark <- mark + 1           
        }
      } 
      if (assoc["muvalue_muvalue",][[m]]) { # muvalue*muvalue
        sel <- assoc["which_interactions",][[m]][["muvalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- invlink_m(eta_m) * invlink_j(eta_j) 
          a_X[[mark]] <- val
          mark <- mark + 1                   
        }
      }       
      # muslope and any interactions
      if (assoc["muslope",][[m]]) { # muslope
        val <- (invlink_m(eps_m) - invlink_m(eta_m)) / eps
        a_X[[mark]] <- val
        mark <- mark + 1                   
      }
      if (assoc["muslope_data",][[m]]) { # muslope*data
        dydt_m <- (invlink_m(eps_m) - invlink_m(eta_m)) / eps
        idx_ms <- which(names(K_data_m) == "muslope_data")
        cbeg  <- sum(K_data_m[0:(idx_ms-1)]) + 1
        cend  <- sum(K_data_m[0: idx_ms   ])
        val   <- as.vector(dydt_m) * data_m[, cbeg:cend, drop = FALSE] 
        a_X[[mark]] <- val
        mark <- mark + 1              
      }    
      # muauc
      if (assoc["muauc",][[m]]) { # muauc
        val   <- c()
        for (j in 1:length(eta_m)) {
          wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
          auc_j <- invlink_m(auc_m[((j-1) * qnodes + 1):(j * qnodes)])
          val[j] <- sum(wgt_j * auc_j)
        }
        a_X[[mark]] <- val
        mark <- mark + 1 
      }
    }
  }
  for (m in 1:M) {
    # shared_b
    if (assoc["shared_b",][[m]]) {
      sel <- assoc["which_b_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1                   
    }
  }    
  for (m in 1:M) {
    # shared_coef
    if (assoc["shared_coef",][[m]]) {
      sel <- assoc["which_coef_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1                   
    }
  }
  return(do.call("cbind", a_X))
}

# Function to get linear predictor and other bits
#
# @param parts A named list containing the parts for constructing the association 
#   structure. It may contain elements $mod_eta, $mod_eps, $mod_auc, etc. as 
#   well as $xmat_data, $K_data
# @param which
get_element <- function(parts, m = 1, which = "eta", ...) {
  dots <- list(...)
  ok_which_args <- c("eta", "eps", "auc", "xmat_data", "K_data", "b_mat")
  if (!which %in% ok_which_args)
    stop("'which' must be one of: ", paste(ok_which_args, collapse = ", "))
  if (which %in% c("eta", "eps", "auc")) {
    part <- parts[[m]][[paste0("mod_", which)]]
    if (is.null(part)) { # model doesn't include an assoc related to 'part'
      return(NULL)
    } else if (!is.null(part$linpred)) { # linpred already provided in object
      return(part$linpred)
    } else { # need to construct linpred
      x <- part$x
      Zt <- part$Zt
      Znames  <- part$Z_names
      if (is.null(x) || is.null(Zt))
        stop(paste0("Bug found: cannot find x and Zt in object. They are ",
             "required to build the linear predictor for '", which, "'."))          
      beta <- dots$beta[[m]]
      b <- dots$b[[m]] 
      b <- pp_b_ord(if (is.matrix(b)) b else t(b), Znames)
      if (is.null(beta) || is.null(b))
        stop("Bug found: beta and b must be provided to construct linpred.")
      return(linear_predictor.default(beta, x) + as.vector(b %*% Zt))
    }
  } else if (which %in% c("xmat_data", "K_data", "b_mat")) {
    return(parts[[m]][[which]])
  } else {
    stop("'which' argument doesn't include a valid entry.")
  }
}

# Function to return the range or SD of the predictors, used for scaling the priors
# This is taken from an anonymous function in stan_glm.fit
#
# @param x A vector
scale_val <- function(x) {
  num.categories <- length(unique(x))
  x.scale <- 1
  if (num.categories == 2) {
    x.scale <- diff(range(x))
  } else if (num.categories > 2) {
    x.scale <- sd(x)
  }
  return(x.scale)
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

# Create "prior.info" attribute needed for prior_summary()
#
# @param user_* The user's priors. These should be passed in after broadcasting 
#   the df/location/scale arguments if necessary.
# @param y_has_intercept Vector of T/F, does each long submodel have an intercept?
# @param y_has_predictors Vector of T/F, does each long submodel have predictors?
# @param y_has_intercept T/F, does event submodel have an intercept?
# @param has_predictors T/F, does event submodel have predictors?
# @param adjusted_prior_*_scale adjusted scales computed if using autoscaled priors
# @param family A list of family objects.
# @return A named list with components 'prior*', 'prior*_intercept', 
#   'prior_covariance' and 'prior*_aux' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_jm_prior <-
  function(user_priorLong = NULL,
           user_priorLong_intercept = NULL,
           user_priorLong_aux = NULL,
           user_priorEvent = NULL,
           user_priorEvent_intercept = NULL,
           user_priorEvent_aux = NULL,
           user_priorAssoc = NULL,
           user_prior_covariance = NULL,
           y_has_intercept = NULL,
           e_has_intercept = NULL,
           y_has_predictors = NULL,
           e_has_predictors = NULL,
           has_assoc = NULL,
           adjusted_priorLong_scale = NULL,
           adjusted_priorLong_intercept_scale = NULL, 
           adjusted_priorLong_aux_scale = NULL,
           adjusted_priorEvent_scale = NULL,
           adjusted_priorEvent_intercept_scale = NULL, 
           adjusted_priorEvent_aux_scale = NULL,           
           adjusted_priorAssoc_scale = NULL,
           family = NULL, 
           basehaz = NULL,
           stub_for_names = "Long") {
    if (!is.null(family) && !is(family, "list"))
      stop("'family' should be a list of family objects, one for each submodel.")
    if (!is.null(has_assoc) && !is.logical(has_assoc) && (length(has_assoc) == 1L))
      stop("'has_assoc' should be a logical vector of length 1.")
    M <- length(family)
    
    prior_list <- list()
    
    if (!is.null(user_priorLong)) {
      rescaled_coefLong <- mapply(check_if_rescaled, user_priorLong, 
                                  y_has_predictors, adjusted_priorLong_scale)
      rescaled_intLong  <- mapply(check_if_rescaled, user_priorLong_intercept, 
                                  y_has_intercept, adjusted_priorLong_intercept_scale)
      rescaled_auxLong  <- mapply(check_if_rescaled, user_priorLong_aux, 
                                  TRUE, adjusted_priorLong_aux_scale) 
      for (m in 1:M) {
        user_priorLong[[m]] <- 
          rename_t_and_cauchy(user_priorLong[[m]], y_has_predictors[m])
        user_priorLong_intercept[[m]] <-
          rename_t_and_cauchy(user_priorLong_intercept[[m]], y_has_intercept[m])
        user_priorLong_aux[[m]] <-
          rename_t_and_cauchy(user_priorLong_aux[[m]], TRUE)
      }
      prior_list$priorLong <- list_nms(lapply(1:M, function(m) {
        if (!y_has_predictors[m]) NULL else with(user_priorLong[[m]], list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefLong[m])
            adjusted_priorLong_scale[[m]] else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))        
      }), M, stub = stub_for_names)
      prior_list$priorLong_intercept <- list_nms(lapply(1:M, function(m) {
        if (!y_has_intercept[m]) NULL else with(user_priorLong_intercept[[m]], list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_intLong[m]) 
            adjusted_priorLong_intercept_scale[[m]] else NULL,
          df = if (prior_dist_name %in% "student_t") 
            prior_df else NULL
        ))
      }), M, stub = stub_for_names)      
      aux_name <- lapply(family, .rename_aux)
      prior_list$priorLong_aux <- list_nms(lapply(1:M, function(m) {
        if (is.na(aux_name[[m]])) NULL else with(user_priorLong_aux[[m]], list(
          dist = prior_dist_name,
          location = if (!is.na(prior_dist_name) && 
                         prior_dist_name != "exponential")
            prior_mean else NULL,
          scale = if (!is.na(prior_dist_name) && 
                      prior_dist_name != "exponential")
            prior_scale else NULL,
          adjusted_scale = if (rescaled_auxLong[m])
            adjusted_priorLong_aux_scale[[m]] else NULL,
          df = if (!is.na(prior_dist_name) && 
                   prior_dist_name %in% "student_t")
            prior_df else NULL, 
          rate = if (!is.na(prior_dist_name) && 
                     prior_dist_name %in% "exponential")
            1 / prior_scale else NULL,
          aux_name = aux_name[[m]]
        ))
      }), M, stub = stub_for_names)     
    }

    if (!is.null(user_priorEvent)) {
      rescaled_coefEvent <- check_if_rescaled(user_priorEvent, e_has_predictors,
                                              adjusted_priorEvent_scale)
      rescaled_intEvent  <- check_if_rescaled(user_priorEvent_intercept, e_has_intercept, 
                                              adjusted_priorEvent_intercept_scale)
      rescaled_auxEvent  <- check_if_rescaled(user_priorEvent_aux, TRUE, 
                                              adjusted_priorEvent_aux_scale)
      user_priorEvent <- 
        rename_t_and_cauchy(user_priorEvent, e_has_predictors)  
      user_priorEvent_intercept <- 
        rename_t_and_cauchy(user_priorEvent_intercept, e_has_intercept)  
      user_priorEvent_aux <- 
        rename_t_and_cauchy(user_priorEvent_aux, TRUE)     
      prior_list$priorEvent <-
        if (!e_has_predictors) NULL else with(user_priorEvent, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefEvent)
            adjusted_priorEvent_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))
      prior_list$priorEvent_intercept <-
        if (!e_has_intercept) NULL else with(user_priorEvent_intercept, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_intEvent)
            adjusted_priorEvent_intercept_scale else NULL,
          df = if (prior_dist_name %in% "student_t")
            prior_df else NULL
        ))
      e_aux_name <- .rename_e_aux(basehaz) 
      prior_list$priorEvent_aux <-
        with(user_priorEvent_aux, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_auxEvent)
            adjusted_priorEvent_aux_scale else NULL,
          df = if (!is.na(prior_dist_name) && 
                   prior_dist_name %in% "student_t")
            prior_df else NULL, 
          aux_name = e_aux_name
        ))      
    }

    if (!is.null(user_priorAssoc)) {
      rescaled_coefAssoc <- check_if_rescaled(user_priorAssoc, has_assoc, 
                                              adjusted_priorAssoc_scale)
      user_priorAssoc <- rename_t_and_cauchy(user_priorAssoc, has_assoc)        
      prior_list$priorAssoc <-
        if (!has_assoc) NULL else with(user_priorAssoc, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefAssoc)
            adjusted_priorAssoc_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))
    }
 
    if (length(user_prior_covariance))
      prior_list$prior_covariance <- user_prior_covariance
    
    return(prior_list)
  }

# Get name of auxiliary parameters for event submodel
#
# @param basehaz A list with information about the baseline hazard
.rename_e_aux <- function(basehaz) {
  nm <- basehaz$type_name
  if (nm == "weibull") "weibull-shape" else
    if (nm == "bs") "spline-coefficients" else
      if (nm == "piecewise") "piecewise-coefficients" else NA
}

# Check if priors were autoscaled
#
# @param prior_stuff A list with prior info returned by handle_glm_prior
# @param has A logical checking, for example, whether the model has_predictors, 
#   has_intercept, has_assoc, etc
# @param adjusted_prior_scale The prior scale after any autoscaling
check_if_rescaled <- function(prior_stuff, has, adjusted_prior_scale) {
  prior_stuff$prior_autoscale && has &&
    !is.na(prior_stuff$prior_dist_name) &&
    !all(prior_stuff$prior_scale == adjusted_prior_scale)      
}

# Rename the t prior as being student-t or cauchy
#
# @param prior_stuff A list with prior info returned by handle_glm_prior
# @param has A logical checking, for example, whether the model has_predictors, 
#   has_intercept, has_assoc, etc
rename_t_and_cauchy <- function(prior_stuff, has) {
  if (has && prior_stuff$prior_dist_name %in% "t") {
    if (all(prior_stuff$prior_df == 1)) {
      prior_stuff$prior_dist_name <- "cauchy"
    } else {
      prior_stuff$prior_dist_name <- "student_t"
    }
  }
  return(prior_stuff)
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
  est <- lapply(y_mod_stuff, function(x) summary(x$mod)$coefficients[, "Estimate"])
  xbar      <- fetch(y_mod_stuff, "xbar")
  gamma     <- lapply(seq_along(y_mod_stuff), function(m) 
    return_intercept(est[[m]]) - xbar[[m]] %*% drop_intercept(est[[m]]))
  gamma_nob <- gamma[as.logical(standata$has_intercept_nob)]
  gamma_lob <- gamma[as.logical(standata$has_intercept_lob)]
  gamma_upb <- gamma[as.logical(standata$has_intercept_upb)]
  beta      <- lapply(est, drop_intercept)
  aux       <- lapply(y_mod_stuff, function(x) sigma(x$mod))
  e_beta    <- e_mod_stuff$mod$coef
  e_aux     <- if (standata$basehaz_type == 1L) runif(1, 0.5, 3) else rep(0, standata$basehaz_df)
  z_beta        <- standardise_coef(unlist(beta), standata$prior_mean,           standata$prior_scale)
  aux_unscaled  <- standardise_coef(unlist(aux),  standata$prior_mean_for_aux,   standata$prior_scale_for_aux)
  aux_unscaled  <- aux_unscaled[as.logical(standata$has_aux)] # only keep aux where relevant
  e_z_beta      <- standardise_coef(e_beta,       standata$e_prior_mean,         standata$e_prior_scale) 
  e_aux_unscaled<- standardise_coef(e_aux,        standata$e_prior_mean_for_aux, standata$e_prior_scale_for_aux)
  b_Cov         <- lapply(y_mod_stuff, function(x) lme4::VarCorr(x$mod)[[1L]])
  sel           <- sapply(y_mod_stuff, function(x) length(lme4::VarCorr(x$mod)) > 1L)
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
  
  # Parameters related to priors
  len_global <- sum((2 * (standata$prior_dist == 3)) + (4 * (standata$prior_dist == 4)))
  len_local2 <- sum((standata$prior_dist == 3) * standata$KM) 
  len_local4 <- sum((standata$prior_dist == 4) * standata$KM)
  len_mix   <- sum((standata$prior_dist %in% c(5,6)) * standata$KM)
  len_ool <- sum(standata$prior_dist == 6)
  len_noise <- sum((standata$family == 8) * standata$NM)
    
  # Function to generate model based initial values
  model_based_inits <- Filter(function(x) (!is.null(x)), list(
    gamma_nob      = array_else_double(gamma_nob),
    gamma_lob      = array_else_double(gamma_lob),
    gamma_upb      = array_else_double(gamma_upb),
    z_beta         = array_else_double(z_beta),
    aux_unscaled   = array_else_double(aux_unscaled),
    e_z_beta       = array_else_double(e_z_beta),
    e_aux_unscaled = array_else_double(e_aux_unscaled),
    e_gamma  = array_else_double(rep(0, standata$e_has_intercept)),
    a_z_beta = array_else_double(rep(0, standata$a_K)),
    z_b      = array_else_double(runif(standata$q, -0.5, 0.5)),
    global   = array_else_double(runif(len_global)),
    e_global = array_else_double(runif(get_nvars_for_hs(standata$e_prior_dist))),
    a_global = array_else_double(runif(get_nvars_for_hs(standata$a_prior_dist))),
    local2   = matrix_of_uniforms(nrow = 2, ncol = len_local2),
    local4   = matrix_of_uniforms(nrow = 4, ncol = len_local4),
    e_local  = matrix_of_uniforms(nrow = get_nvars_for_hs(standata$e_prior_dist), ncol = standata$e_K),
    a_local  = matrix_of_uniforms(nrow = get_nvars_for_hs(standata$a_prior_dist), ncol = standata$a_K),
    mix   = if (len_mix > 0) matrix(rep(1, len_mix), 1, len_mix) else matrix(0,0,0),
    e_mix = if (standata$e_prior_dist %in% c(5,6)) matrix(rep(1, standata$e_K), 1, standata$e_K) else matrix(0,0,standata$e_K),
    a_mix = if (standata$a_prior_dist %in% c(5,6)) matrix(rep(1, standata$a_K), 1, standata$a_K) else matrix(0,0,standata$a_K),
    ool   = if (len_ool > 0) as.array(len_ool) else as.array(double(0)), 
    e_ool = if (standata$e_prior_dist == 6) as.array(1) else as.array(double(0)), 
    a_ool = if (standata$a_prior_dist == 6) as.array(1) else as.array(double(0)),
    noise = if (len_noise > 0) matrix(runif(len_noise), 1, len_noise) else matrix(0,0,0),
    z_T   = array_else_double(rep(sqrt(1 / len_z_T), len_z_T)),
    rho   = array_else_double(rep(1 / (len_rho + 1), len_rho)),
    zeta  = array_else_double(normalised_zetas),
    tau   = array_else_double(tau)))
  
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
  default_max_treedepth <- 11L
  
  if (!is.null(user_adapt_delta))
    args$control$adapt_delta <- user_adapt_delta else 
      if (is.null(args$control$adapt_delta))
        args$control$adapt_delta <- default_adapt_delta
  
  if (!is.null(user_max_treedepth))
    args$control$max_treedepth <- user_max_treedepth else
      if (is.null(args$control$max_treedepth))
        args$control$max_treedepth <- default_max_treedepth
  
  if (!"iter" %in% unms) args$iter <- 1000
  if (!"chains" %in% unms) args$chains <- 3
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
    if (!all(check)) STOP_arg(deparse(substitute(arg)), type, null_ok = null_ok)
  } else {
    STOP_arg(deparse(substitute(arg)), type, null_ok = null_ok)
  }
  if (!is.null(validate_length)) {
    if (length(arg) == 1L && validate_length > 1) {
      arg <- if (broadcast) rep(arg, times = validate_length) else arg
    } else if (!length(arg) == validate_length) {
      stop(paste(deparse(substitute(arg)), "is a list of the incorrect length."),
           call. = FALSE)
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
         "longitudinal submodels, but the list is of the incorrect length.", 
         call. = FALSE)
  }
}

# From a vector of length M giving the number of elements (for example number
# of parameters or observations) for each submodel, create an indexing array 
# of dimension M * 2, where column 1 is the beginning index and 2 is the end index
#
# @param x A numeric vector
# @return A length(x) * 2 array
get_idx_array <- function(x) {
  as.array(do.call("rbind", lapply(1:length(x), function(i) {
    idx_beg <- ifelse(x[i] > 0L, sum(x[0:(i-1)]) + 1, 0L)
    idx_end <- ifelse(x[i] > 0L, sum(x[0:i]),         0L)
    c(idx_beg, idx_end)
  })))
}

# Error message when the argument contains an object of the incorrect type
STOP_arg <- function(arg_name, type, null_ok = FALSE) {
  stop(paste0("'", arg_name, "' should be ", ifelse(null_ok, "NULL, ", ""), "a ", 
              type, " object or a list of ", type, " objects."), call. = FALSE) 
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
# Wrapper for using fetch with unlist = TRUE and 
# returning array. Also converts logical to integer.
fetch_array <- function(x, y) {
  val <- fetch(x, y, unlist = TRUE)
  if (is.logical(val)) val <- as.integer(val)
  as.array(val)
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




