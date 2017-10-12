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
#' @param grp_assoc Character string specifying the method for combining information
#'   across the lower level units clustered within an individual when forming the
#'   association structure. This is only relevant when a grouping factor is  
#'   specified in \code{formulaLong} that corresponds to clustering within 
#'   individuals. Can be \code{"sum"} for specifying which indicates
#'   the association structure should be based on a summation across the lower
#'   level units clustered within an individual. Can be \code{"mean"} which
#'   indicates that the association structure should be based on the mean
#'   (i.e. average) taken across the lower level units clustered within an
#'   individual. (As an example, specifying \code{assoc = "muvalue"} 
#'   and \code{grp_assoc = "sum"} would mean the log hazard 
#'   for an individual would be linearly related to the sum of
#'   the expected values for each of the lower level units clustered within
#'   that individual). 
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
#' @param epsilon The half-width of the central difference used to numerically
#'   calculate the derivate when the \code{"etaslope"} association structure 
#'   is used.   
#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
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
#'   number of quadrature nodes, specified through the \code{qnodes} 
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
#'               formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
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
#'               formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
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
#'               formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
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
#'               formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
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
#'         formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
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
#'         formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
#'         dataEvent = pbcSurv,
#'         assoc = list(c("etavalue", "etavalue_etavalue(2)"), "etavalue"),
#'         time_var = "year", 
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' }
#'  
#' @importFrom lme4 lmerControl glmerControl glmer
#' 
stan_jm <- function(formulaLong, dataLong, formulaEvent, dataEvent, time_var, 
                    id_var, family = gaussian, assoc = "etavalue", 
                    lag_assoc = 0, grp_assoc, dataAssoc, epsilon = 1E-5,
                    basehaz = c("weibull", "bs", "piecewise"), basehaz_ops, 
                    qnodes = 15, init = "prefit", weights, ...,	
                    priorLong = normal(), priorLong_intercept = normal(), 
                    priorLong_aux = cauchy(0, 5), priorEvent = normal(), 
                    priorEvent_intercept = normal(), priorEvent_aux = cauchy(),
                    priorEvent_assoc = normal(), prior_covariance = decov(), 
                    prior_PD = FALSE, algorithm = c("sampling", "meanfield", "fullrank"), 
                    adapt_delta = NULL, max_treedepth = NULL, QR = FALSE, 
                    sparse = FALSE) {
  
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------  
  
  # Set seed if specified
  dots <- list(...)
  if ("seed" %in% names(dots))
    set.seed(dots$seed)
  
  algorithm <- match.arg(algorithm)
  basehaz   <- match.arg(basehaz)
  
  if (missing(offset))      offset      <- NULL 
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  if (missing(weights))     weights     <- NULL
  if (missing(id_var))      id_var      <- NULL
  if (missing(time_var))    time_var    <- NULL
  if (missing(grp_assoc))   grp_assoc   <- NULL
  if (missing(dataAssoc))   dataAssoc   <- NULL 
  
  if (!is.null(weights)) 
    stop("'weights' are not yet implemented.")
  if (!is.null(dataAssoc))
    stop("'dataAssoc' argument not yet implemented.")
  if (QR)               
    stop("'QR' decomposition is not yet implemented.")
  if (sparse)
    stop("'sparse' option is not yet implemented.")
  
  if (is.null(time_var))
    stop("'time_var' must be specified.")
  if (is.null(id_var))
    stop("'id_var' must be specified.")

  # Formula
  formulaLong <- validate_arg(formulaLong, "formula"); M <- length(formulaLong)
  
  # Data
  dataLong <- validate_arg(dataLong, "data.frame", validate_length = M)  
  dataLong <- check_vars_are_included(dataLong, id_var, time_var)
  dataLong <- xapply(formulaLong, dataLong, FUN = get_all_vars)
  
  dataEvent <- as.data.frame(dataEvent)
  dataEvent <- check_vars_are_included(dataEvent, id_var)
  dataEvent <- get_all_vars(formulaEvent, dataEvent, dataEvent[[id_var]])
  names(dataEvent) <- c(names(dataEvent), id_var)
    
  # Family
  ok_family_classes <- c("function", "family", "character")
  ok_families <- c("binomial", "gaussian", "Gamma", 
                   "inverse.gaussian", "poisson", "neg_binomial_2")
  family <- validate_arg(family, ok_family_classes, validate_length = M)
  family <- lapply(family, validate_famlink, ok_families)
  
  # Assoc
  ok_assoc_classes <- c("NULL", "character")
  assoc <- validate_arg(assoc, ok_assoc_classes, validate_length = M)

  # Is priorLong* already a list?
  priorLong <- broadcast_prior(priorLong, M)
  priorLong_intercept <- broadcast_prior(priorLong_intercept, M)
  priorLong_aux <- broadcast_prior(priorLong_aux, M)
 
  # Observation weights
  if (!is.null(weights)) 
    weights <- check_weights(weights, id_var)
  
  #-----------
  # Fit model
  #----------- 
  
  stanfit <- stan_jm.fit(formulaLong = formulaLong, dataLong = dataLong, 
                         formulaEvent = formulaEvent, dataEvent = dataEvent, 
                         time_var = time_var, id_var = id_var, family = family,
                         assoc = assoc, lag_assoc = lag_assoc, grp_assoc = grp_assoc, 
                         dataAssoc = dataAssoc, epsilon = epsilon, basehaz = basehaz, 
                         basehaz_ops = basehaz_ops, qnodes = qnodes, init = init, 
                         weights = weights, ..., priorLong = priorLong, 
                         priorLong_intercept = priorLong_intercept, 
                         priorLong_aux = priorLong_aux, priorEvent = priorEvent, 
                         priorEvent_intercept = priorEvent_intercept, 
                         priorEvent_aux = priorEvent_aux, priorEvent_assoc = priorEvent_assoc, 
                         prior_covariance = prior_covariance, prior_PD = prior_PD, 
                         algorithm = algorithm, adapt_delta = adapt_delta, 
                         max_treedepth = max_treedepth, QR = QR, sparse = sparse)

  y_mod <- attr(stanfit, "y_mod")
  e_mod <- attr(stanfit, "e_mod")
  cnms  <- attr(stanfit, "cnms")
  basehaz <- attr(stanfit, "basehaz")
  stanfit <- drop_attributes(stanfit, "y_mod", "e_mod", "cnms", "basehaz")
  
  terms <- c(fetch(y_mod, "terms"), list(e_mod$terms))
  n_yobs <- fetch_(y_mod, "X", "N")
  n_pats <- e_mod$Npat
  #n_grps <- standata$l - 1
  #names(n_grps) <- cnms_nms  # n_grps is num. of levels within each grouping factor
  #names(p) <- cnms_nms       # p is num. of variables within each grouping factor

  fit <- nlist(stanfit, formula = c(formulaLong, formulaEvent), family, terms,
               id_var, time_var, weights, qnodes, basehaz,
               M, cnms, n_grps, n_pats, n_yobs, 
               
               assoc, a_mod_stuff, clust_stuff,
               grp_assoc = if (any(unlist(has_clust))) grp_assoc else NULL,
               fr = list_nms(c(fetch(a_mod_stuff, "model_frame"), 
                               list(e_mod_stuff$model_frame)), M),
               
               algorithm, terms, glmod = y_mod, survmod = e_mod,
               model_dataLong = dataLong, model_dataEvent = dataEvent,
               prior.info = NULL, stan_function = "stan_jm", 
               call = match.call(expand.dots = TRUE))
  
  out <- stanmvreg(fit)
  return(out)
}

  #--------------------------------
  # Data for association structure
  #--------------------------------
  
 

  # Return design matrices for evaluating longitudinal submodel quantities
  # at the quadrature points
  eps <- 1E-5 # time shift for numerically calculating deriv using one-sided diff
  auc_qnodes <- 15L
  a_mod_stuff <- mapply(handle_assocmod, 1:M, m_mc, dataLong, y_mod_stuff,
                        clust_stuff = clust_stuff, SIMPLIFY = FALSE, 
                        MoreArgs = list(id_list         = e_mod_stuff$cids, 
                                        times           = e_mod_stuff$cpts, 
                                        assoc           = assoc, 
                                        id_var          = id_var, 
                                        time_var        = time_var, 
                                        eps             = eps, 
                                        auc_qnodes      = auc_qnodes,
                                        dataAssoc       = dataAssoc,
                                        env             = calling_env))

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
    a_xbar                     = if (a_K) a_prior_stuff$a_xbar else numeric(0),
	
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
  standata$Npat      <- as.integer(e_mod_stuff$Npat)
  standata$Nevents   <- as.integer(e_mod_stuff$Nevents)
  standata$qnodes    <- as.integer(qnodes)
  standata$qwts      <- as.array(e_mod_stuff$qwts)
  standata$Npat_times_qnodes <- as.integer(e_mod_stuff$Npat * qnodes)
  standata$e_times <- c(e_mod_stuff$eventtime[e_mod_stuff$status == 1], unlist(e_mod_stuff$qpts))
  standata$nrow_e_Xq <- length(standata$e_times)
  standata$e_has_intercept <- as.integer(e_has_intercept)
  standata$e_Xq      <- e_mod_stuff$xtemp
  standata$e_xbar    <- as.array(e_mod_stuff$xbar)
  standata$e_weights       <- as.array(e_weights)
  standata$e_weights_rep   <- as.array(rep(e_weights, times = qnodes))
  
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
  standata$norm_const <- e_mod_stuff$norm_const
  
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
  
  # Data for association structure when there is
  # clustering below the patient-level
  standata$has_clust <- as.array(as.integer(has_clust))
  if (any(has_clust)) { # has lower level clustering
    parts_clust_mat <- rstan::extract_sparse_parts(clust_mat[[1L]])
    standata$clust_nnz <- length(parts_clust_mat$w)
    standata$clust_w <- parts_clust_mat$w
    standata$clust_v <- parts_clust_mat$v
    standata$clust_u <- parts_clust_mat$u
  } else { # no lower level clustering
    standata$clust_nnz <- 0L
    standata$clust_w <- double(0)
    standata$clust_v <- integer(0)
    standata$clust_u <- integer(0)
  }
  
  # Data for calculating value, slope, auc in GK quadrature 
  standata$nrow_y_Xq <- as.array(as.integer(
    sapply(a_mod_stuff, function(x) NROW(x$mod_eta$xtemp))))
  standata$idx_q <- get_idx_array(standata$nrow_y_Xq)
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
  standata$auc_qnodes <- 
    as.integer(auc_qnodes)
  standata$Npat_times_auc_qnodes <- 
    as.integer(e_mod_stuff$Npat * auc_qnodes) 
  standata$nrow_y_Xq_auc <- as.array(as.integer(
    sapply(a_mod_stuff, function(x) NROW(x$mod_auc$xtemp))))
  standata$idx_qauc <- get_idx_array(standata$nrow_y_Xq_auc)
  auc_qwts <- 
    unlist(lapply(e_mod_stuff$qtimes, function(x) 
      lapply(x, function(y) 
        lapply(get_quadpoints(auc_qnodes)$weights, unstandardise_qwts, 0, y)))) 
  standata$auc_qwts <- 
    if (standata$assoc_uses[3]) as.array(auc_qwts) else double(0)

  # Interactions between association terms and data
  # design matrix for the interactions
  standata$y_Xq_data <- 
    as.array(as.matrix(Matrix::bdiag(fetch(a_mod_stuff, "xmat_data"))))
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
    dropargs <- c("chains", "cores", "iter", "refresh")
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
  





