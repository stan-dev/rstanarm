# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2018 Sam Brilleman
# Copyright (C) 2018 Trustees of Columbia University
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

#' Bayesian survival models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for survival models (sometimes known as models for 
#' time-to-event data). Currently, the command fits:
#' (i) flexible parametric (cubic spline-based) survival 
#' models on the hazard scale, with covariates included under assumptions of 
#' either proportional or non-proportional hazards;
#' (ii) standard parametric (exponential, Weibull and Gompertz) survival 
#' models on the hazard scale, with covariates included under assumptions of 
#' either proportional or non-proportional hazards; and
#' (iii) standard parametric (exponential, Weibull) accelerated failure time
#' models, with covariates included under assumptions of either time-fixed or 
#' time-varying survival time ratios. Left, right, and interval censored 
#' survival data are allowed. Delayed entry is allowed. Both fixed and random
#' effects can be estimated for covariates (i.e. group-specific parameters
#' are allowed). Time-varying covariates and time-varying coefficients are 
#' both allowed. For modelling each time-varying coefficient (i.e. time-varying 
#' log hazard ratio or time-varying log survival time ratio) the user can 
#' choose between either a smooth B-spline function or a piecewise constant  
#' function.
#'
#' @export
#' @importFrom splines bs
#' @import splines2
#'  
#' @template args-dots
#' @template args-priors
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' 
#' @param formula A two-sided formula object describing the model. 
#'   The left hand side of the formula should be a \code{Surv()} 
#'   object. Left censored, right censored, and interval censored data 
#'   are allowed, as well as delayed entry (i.e. left truncation). See 
#'   \code{\link[survival]{Surv}} for how to specify these outcome types. 
#'   The right hand side of the formula can include fixed and/or random
#'   effects of covariates, with random effects specified in the same 
#'   way as for the \code{\link[lme4]{lmer}} function in the \pkg{lme4}
#'   package. If you wish to include time-varying effects (i.e. time-varying 
#'   coefficients, e.g. non-proportional hazards) in the model
#'   then any covariate(s) that you wish to estimate a time-varying 
#'   coefficient for should be specified as \code{tve(varname)} where 
#'   \code{varname} is the name of the covariate. For more information on 
#'   how time-varying effects are formulated see the documentation
#'   for the \code{\link{tve}} function as well as the \strong{Details} 
#'   and \strong{Examples} sections below.
#' @param data A data frame containing the variables specified in 
#'   \code{formula}.
#' @param basehaz A character string indicating which baseline hazard or
#'   baseline survival distribution to use for the event submodel. 
#'   
#'   The following are available under a hazard scale formulation: 
#'   \itemize{
#'     \item \code{"ms"}: A flexible parametric model using cubic M-splines to 
#'     model the baseline hazard. The default locations for the internal knots, 
#'     as well as the basis terms for the splines, are calculated with respect
#'     to time. If the model does \emph{not} include any time-dependendent 
#'     effects then a closed form solution is available for both the hazard
#'     and cumulative hazard and so this approach should be relatively fast.
#'     On the other hand, if the model does include time-varying effects then
#'     quadrature is used to evaluate the cumulative hazard at each MCMC
#'     iteration and, therefore, estimation of the model will be slower.
#'     \item \code{"bs"}: A flexible parametric model using cubic B-splines to 
#'     model the \emph{log} baseline hazard. The default locations for the  
#'     internal knots, as well as the basis terms for the splines, are calculated 
#'     with respect to time. A closed form solution for the cumulative hazard 
#'     is \strong{not} available regardless of whether or not the model includes
#'     time-varying effects; instead, quadrature is used to evaluate 
#'     the cumulative hazard at each MCMC iteration. Therefore, if your model
#'     does not include any time-varying effects, then estimation using the 
#'     \code{"ms"} baseline hazard will be faster.
#'     \item \code{"exp"}: An exponential distribution for the event times
#'     (i.e. a constant baseline hazard).
#'     \item \code{"weibull"}: A Weibull distribution for the event times.
#'     \item \code{"gompertz"}: A Gompertz distribution for the event times.
#'   }
#'   
#'   The following are available under an accelerated failure time (AFT)
#'   formulation: 
#'   \itemize{
#'     \item \code{"exp-aft"}: an exponential distribution for the event times.
#'     \item \code{"weibull-aft"}: a Weibull distribution for the event times.
#'   }
#' @param basehaz_ops A named list specifying options related to the baseline
#'   hazard. Currently this can include: \cr
#'   \itemize{
#'     \item \code{degree}: A positive integer specifying the degree for the 
#'     M-splines or B-splines. The default is \code{degree = 3}, which
#'     corresponds to cubic splines. Note that specifying \code{degree = 0}
#'     is also allowed and corresponds to piecewise constant.
#'     \item \code{df}: A positive integer specifying the degrees of freedom 
#'     for the M-splines or B-splines. For M-splines (i.e. when 
#'     \code{basehaz = "ms"}), two boundary knots and \code{df - degree - 1} 
#'     internal knots are used to generate the spline basis. For B-splines 
#'     (i.e. when \code{basehaz = "bs"}), two boundary knots and 
#'     \code{df - degree} internal knots are used to generate the spline 
#'     basis. The difference is due to the fact that the M-spline basis
#'     includes an intercept, whereas the B-spline basis does not. The 
#'     default is \code{df = 6} for M-splines and \code{df = 5} for
#'     B-splines (i.e. two boundary knots and two internal knots when the
#'     default cubic splines are being used). The internal knots are placed 
#'     at equally spaced percentiles of the distribution of uncensored event 
#'     times.
#'     \item \code{knots}: A numeric vector explicitly specifying internal 
#'     knot locations for the M-splines or B-splines. Note that \code{knots} 
#'     cannot be specified if \code{df} is specified.
#'   }
#'   Note that for the M-splines and B-splines -- in addition to any internal
#'   \code{knots} -- a lower boundary knot is placed at the earliest entry time
#'   and an upper boundary knot is placed at the latest event or censoring time.
#'   These boundary knot locations are the default and cannot be changed by the
#'   user.
#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard when \code{basehaz = "bs"}
#'   or when time-varying effects are specified in the linear predictor. 
#'   Options are 15 (the default), 11 or 7.
#' @param prior_intercept The prior distribution for the intercept in the 
#'   linear predictor. All models include an intercept parameter.
#'   \code{prior_intercept} can be a call to \code{normal}, 
#'   \code{student_t} or \code{cauchy}. See the \link[=priors]{priors help page} 
#'   for details on these functions. However, note that default scale for 
#'   \code{prior_intercept} is 20 for \code{stan_surv} models (rather than 10,
#'   which is the default scale used for \code{prior_intercept} by most 
#'   \pkg{rstanarm} modelling functions). To omit a prior on the intercept 
#'   ---i.e., to use a flat (improper) uniform prior--- \code{prior_intercept} 
#'   can be set to \code{NULL}.
#'   
#'   \strong{Note:} The prior distribution for the intercept is set so it
#'   applies to the value \emph{when all predictors are centered} and with an  
#'   adjustment ("constant shift") equal to the \emph{log crude event rate}.
#'   However, the reported \emph{estimates} for the intercept always correspond 
#'   to a parameterization without centered predictors and without the 
#'   "constant shift". That is, these adjustments are made internally to help
#'   with numerical stability and sampling, but the necessary 
#'   back-transformations are made so that they are not relevant for the 
#'   estimates returned to the user.
#' @param prior_aux The prior distribution for "auxiliary" parameters related 
#'   to the baseline hazard. The relevant parameters differ depending 
#'   on the type of baseline hazard specified in the \code{basehaz} 
#'   argument. The following applies (for further technical details, 
#'   refer to the \emph{stan_surv: Survival (Time-to-Event) Models vignette}):
#'   \itemize{
#'     \item \code{basehaz = "ms"}: the auxiliary parameters are the 
#'     coefficients for the M-spline basis terms on the baseline hazard. 
#'     These coefficients are defined using a simplex; that is, they are 
#'     all between 0 and 1, and constrained to sum to 1. This constraint 
#'     is necessary for identifiability of the intercept in the linear 
#'     predictor. The default prior is a Dirichlet distribution with all 
#'     concentration parameters set equal to 1. That is, a uniform 
#'     prior over all points defined within the support of the simplex. 
#'     Specifying all concentration parameters equal and > 1 supports a more 
#'     even distribution (i.e. a smoother spline function), while specifying a 
#'     all concentration parameters equal and < 1 supports a more sparse 
#'     distribution (i.e. a less smooth spline function).
#'     \item \code{basehaz = "bs"}: the auxiliary parameters are the 
#'     coefficients for the B-spline basis terms on the log baseline hazard. 
#'     These parameters are unbounded. The default prior is a normal 
#'     distribution with mean 0 and scale 20.
#'     \item \code{basehaz = "exp"} or \code{basehaz = "exp-aft"}: 
#'     there is \strong{no} auxiliary parameter,
#'     since the log scale parameter for the exponential distribution is 
#'     incorporated as an intercept in the linear predictor.
#'     \item \code{basehaz = "weibull"} or \code{basehaz = "weibull-aft"}: 
#'     the auxiliary parameter is the Weibull 
#'     shape parameter, while the log scale parameter for the Weibull 
#'     distribution is incorporated as an intercept in the linear predictor.
#'     The auxiliary parameter has a lower bound at zero. The default prior is  
#'     a half-normal distribution with mean 0 and scale 2.
#'     \item \code{basehaz = "gompertz"}: the auxiliary parameter is the Gompertz 
#'     scale parameter, while the log shape parameter for the Gompertz 
#'     distribution is incorporated as an intercept in the linear predictor.
#'     The auxiliary parameter has a lower bound at zero. The default prior is  
#'     a half-normal distribution with mean 0 and scale 0.5.
#'   }
#'   Currently, \code{prior_aux} can be a call to \code{dirichlet}, 
#'   \code{normal}, \code{student_t}, \code{cauchy} or \code{exponential}. 
#'   See \code{\link{priors}} for details on these functions. Note that not 
#'   all prior distributions are allowed with all types of baseline hazard. 
#'   To omit a prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_aux} to \code{NULL}. 
#' @param prior_smooth This is only relevant when time-varying effects are 
#'   specified in the model (i.e. the \code{tve()} function is used in the 
#'   model formula. When that is the case, \code{prior_smooth} determines the
#'   prior distribution given to the hyperparameter (standard deviation) 
#'   contained in a random-walk prior for the parameters of the function 
#'   used to generate the time-varying coefficient (i.e. the B-spline
#'   coefficients when a B-spline function is used to model the time-varying
#'   coefficient, or the deviations in the log hazard ratio specific to each
#'   time interval when a piecewise constant function is used to model the 
#'   time-varying coefficient). Lower values for the hyperparameter
#'   yield a less flexible function for the time-varying coefficient. 
#'   Specifically, \code{prior_smooth} can be a call to \code{exponential} to 
#'   use an exponential distribution, or \code{normal}, \code{student_t} or 
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy 
#'   prior. See \code{\link{priors}} for details on these functions. To omit a 
#'   prior ---i.e., to use a flat (improper) uniform prior--- set 
#'   \code{prior_smooth} to \code{NULL}. The number of hyperparameters depends
#'   on the model specification (i.e. the number of time-varying effects
#'   specified in the model) but a scalar prior will be recycled as necessary
#'   to the appropriate length.
#'  
#' @details
#' \subsection{Model formulations}{
#'   Let \eqn{h_i(t)} denote the hazard for individual \eqn{i} at time 
#'   \eqn{t}, \eqn{h_0(t)} the baseline hazard at time \eqn{t}, \eqn{X_i} 
#'   a vector of covariates for individual \eqn{i}, \eqn{\beta} a vector of 
#'   coefficients, \eqn{S_i(t)} the survival probability for individual 
#'   \eqn{i} at time \eqn{t}, and \eqn{S_0(t)} the baseline survival 
#'   probability at time \eqn{t}. Without time-varying effects in the 
#'   model formula our linear predictor is \eqn{\eta_i = X_i \beta}, whereas
#'   with time-varying effects in the model formula our linear predictor
#'   is \eqn{\eta_i(t) = X_i(t) \beta(t)}. Then the following definitions of 
#'   the hazard function and survival function apply:
#'   
#'   \tabular{llll}{
#'     \strong{Scale    }                                                \tab 
#'     \strong{TVE      }                                                \tab
#'     \strong{Hazard   }                                                \tab 
#'     \strong{Survival }                                                \cr
#'     \emph{Hazard}                                                     \tab 
#'     \emph{No}                                                         \tab
#'     \eqn{h_i(t) = h_0(t) \exp(\eta_i)}                                \tab
#'     \eqn{S_i(t) = [S_0(t)]^{\exp(\eta_i)}}                            \cr 
#'     \emph{Hazard}                                                     \tab 
#'     \emph{Yes}                                                        \tab
#'     \eqn{h_i(t) = h_0(t) \exp(\eta_i(t))}                             \tab
#'     \eqn{S_i(t) = \exp(- \int_0^t h_i(u) du )}                        \cr
#'     \emph{AFT}                                                        \tab 
#'     \emph{No}                                                         \tab
#'     \eqn{h_i(t) = \exp(-\eta_i) h_0 (t \exp(-\eta_i))}                \tab
#'     \eqn{S_i(t) = S_0 ( t \exp(-\eta_i) )}                            \cr     
#'     \emph{AFT}                                                        \tab 
#'     \emph{Yes}                                                        \tab
#'     \eqn{h_i(t) = \exp(-\eta_i(t)) h_0(\int_0^t \exp(-\eta_i(u)) du)} \tab
#'     \eqn{S_i(t) = S_0 (\int_0^t \exp(-\eta_i(u)) du)}                 \cr
#'   }
#'   
#'   where \emph{AFT} stands for an accelerated failure time formulation, 
#'   and \emph{TVE} stands for time-varying effects in the model formula.
#'   
#'   For models without time-varying effects, the value of \eqn{S_i(t)} can
#'   be calculated analytically (with the one exception being when B-splines 
#'   are used to model the log baseline hazard, i.e. \code{basehaz = "bs"}).
#'   
#'   For models with time-varying effects \eqn{S_i(t)} cannot be calculated 
#'   analytically and so Gauss-Kronrod quadrature is used to approximate the 
#'   relevant integral. The number of nodes used in the quadrature can be 
#'   controlled via the \code{nodes} argument.
#'   
#'   For models estimated on the hazard scale, a hazard ratio can be calculated 
#'   as \eqn{\exp(\beta)}. For models estimated on the AFT scale, a survival 
#'   time ratio can be calculated as \eqn{\exp(\beta)} and an acceleration 
#'   factor can be calculated as \eqn{\exp(-\beta)}.
#'   
#'   Note that the \emph{stan_surv: Survival (Time-to-Event) Models} vignette 
#'   provides more extensive details on the model formulations, including the
#'   parameterisations for each of the parametric distributions.
#' }
#' \subsection{Time-varying effects}{
#'   By default, any covariate effects specified in the \code{formula} are
#'   included in the model under a proportional hazards assumption (for models
#'   estimated using a hazard scale formulation) or under the assumption of
#'   time-fixed acceleration factors (for models estimated using an accelerated
#'   failure time formulation).
#'   
#'   To relax this assumption, it is possible to 
#'   estimate a time-varying effect (i.e. a time-varying coefficient) for a 
#'   given covariate. A time-varying effect is specified in the model 
#'   \code{formula} by wrapping the covariate name in the \code{\link{tve}} 
#'   function. 
#'   
#'   The following applies:
#'   
#'   \itemize{
#'   \item Estimating a time-varying effect within a hazard scale model 
#'   formulation (i.e. when \code{basehaz} is set equal to \code{"ms"}, 
#'   \code{"bs"}, \code{"exp"}, \code{"weibull"} or \code{"gompertz"}) leads
#'   to the estimation of a time-varying hazard ratio for the relevant 
#'   covariate (i.e. non-proportional hazards).
#'   \item Estimating a time-varying effect within an accelerated failure 
#'   time model formulation (i.e. when \code{basehaz} is set equal to 
#'   \code{"exp-aft"}, or \code{"weibull-aft"}) leads to the estimation of a 
#'   time-varying survival time ratio -- or equivalently, a time-varying 
#'   acceleration factor -- for the relevant covariate.
#'   }
#'   
#'   For example, if we wish to estimate a time-varying effect for the 
#'   covariate \code{sex} then we can specify \code{tve(sex)} in the 
#'   \code{formula}, e.g. \code{Surv(time, status) ~ tve(sex) + age + trt}. 
#'   The coefficient for \code{sex} will then be modelled using a flexible 
#'   smooth function based on a cubic B-spline expansion of time.
#'   
#'   Additional arguments used to control the modelling of the time-varying 
#'   effect are explained in the \code{\link{tve}} documentation.
#'   Of particular note is the fact that a piecewise constant basis is 
#'   allowed as a special case of the B-splines. For example, specifying
#'   \code{tve(sex, degree = 0)} in the model formula instead of just
#'   \code{tve(sex)} would request a piecewise constant time-varying effect.
#'   The user can also control the degrees of freedom or knot locations for
#'   the B-spline (or piecewise constant) function.
#'   
#'   It is worth noting that an additional way to control the
#'   flexibility of the function used to model the time-varying effect
#'   is through priors. A random walk prior is used for the piecewise 
#'   constant or B-spline coefficients, and the hyperparameter (standard 
#'   deviation) of the random walk prior can be controlled via the 
#'   \code{prior_smooth} argument. This is a more indirect way to 
#'   control the "smoothness" of the function used to model the time-varying
#'   effect, but it nonetheless might be useful in some settings. The
#'   \emph{stan_surv: Survival (Time-to-Event) Models} vignette provides
#'   more explicit details on the formulation of the time-varying effects
#'   and the prior distributions used for their coefficients.
#'   
#'   It is worth noting that reliable estimation of a time-varying effect  
#'   usually requires a relatively large number of events in the data (e.g. 
#'   say >1000 depending on the setting).
#' }
#'              
#' @examples
#' \donttest{
#' #----- Proportional hazards
#' 
#' # Simulated data
#' library(simsurv)
#' covs <- data.frame(id  = 1:200, 
#'                    trt = stats::rbinom(200, 1L, 0.5))
#' d1 <- simsurv(lambdas = 0.1, 
#'               gammas  = 1.5, 
#'               betas   = c(trt = -0.5),
#'               x       = covs, 
#'               maxt    = 5)
#' d1 <- merge(d1, covs)
#' f1 <- Surv(eventtime, status) ~ trt
#' m1a <- stan_surv(f1, d1, basehaz = "ms",       chains=1,refresh=0,iter=600)
#' m1b <- stan_surv(f1, d1, basehaz = "exp",      chains=1,refresh=0,iter=600)
#' m1c <- stan_surv(f1, d1, basehaz = "weibull",  chains=1,refresh=0,iter=600)
#' m1d <- stan_surv(f1, d1, basehaz = "gompertz", chains=1,refresh=0,iter=600)
#' get_est <- function(x) { fixef(x)["trt"] }
#' do.call(rbind, lapply(list(m1a, m1b, m1c, m1d), get_est))
#' bayesplot::bayesplot_grid(plot(m1a), # compare baseline hazards 
#'                           plot(m1b), 
#'                           plot(m1c), 
#'                           plot(m1d), 
#'                           ylim = c(0, 0.8))
#'     
#' #----- Left and right censored data
#' 
#' # Mice tumor data
#' m2 <- stan_surv(Surv(l, u, type = "interval2") ~ grp, 
#'                 data = mice, chains = 1, refresh = 0, iter = 600)
#' print(m2, 4)
#' 
#' #----- Non-proportional hazards - B-spline tve()
#' 
#' # Simulated data
#' library(simsurv)
#' covs <- data.frame(id  = 1:250, 
#'                    trt = stats::rbinom(250, 1L, 0.5))
#' d3 <- simsurv(lambdas = 0.1, 
#'               gammas  = 1.5, 
#'               betas   = c(trt = -0.5),
#'               tve     = c(trt = 0.2),
#'               x       = covs, 
#'               maxt    = 5)
#' d3 <- merge(d3, covs)
#' m3 <- stan_surv(Surv(eventtime, status) ~ tve(trt), 
#'                 data = d3, chains = 1, refresh = 0, iter = 600)
#' print(m3, 4)
#' plot(m3, "tve") # time-varying hazard ratio
#' 
#' #----- Non-proportional hazards - piecewise constant tve()
#' 
#' # Simulated data
#' library(simsurv)
#' covs <- data.frame(id  = 1:250, 
#'                    trt = stats::rbinom(250, 1L, 0.5))
#' d4 <- simsurv(lambdas = 0.1, 
#'               gammas  = 1.5, 
#'               betas   = c(trt = -0.5),
#'               tve     = c(trt = 0.4),
#'               tvefun  = function(t) { (t > 2.5) },
#'               x       = covs, 
#'               maxt    = 5)
#' d4 <- merge(d4, covs)
#' m4 <- stan_surv(Surv(eventtime, status) ~ 
#'                   tve(trt, degree = 0, knots = c(2.5)), 
#'                 data = d4, chains = 1, refresh = 0, iter = 600)
#' print(m4, 4)
#' plot(m4, "tve") # time-varying hazard ratio
#' 
#' #---------- Compare PH and AFT parameterisations
#' 
#' # Breast cancer data
#' sel <- sample(1:nrow(bcancer), 100)
#' 
#' m_ph  <- stan_surv(Surv(recyrs, status) ~ group, 
#'                    data    = bcancer[sel,], 
#'                    basehaz = "weibull", 
#'                    chains  = 1,
#'                    refresh = 0,
#'                    iter    = 600,
#'                    seed    = 123)
#' m_aft <- stan_surv(Surv(recyrs, status) ~ group, 
#'                    data    = bcancer[sel,], 
#'                    basehaz = "weibull-aft", 
#'                    chains  = 1,
#'                    refresh = 0,
#'                    iter    = 600,
#'                    seed    = 123)
#'
#' exp(fixef(m_ph)) [c('groupMedium', 'groupPoor')] # hazard ratios
#' exp(fixef(m_aft))[c('groupMedium', 'groupPoor')] # survival time ratios
#' 
#' # same model (...slight differences due to sampling)
#' summary(m_ph,  par = "log-posterior")[, 'mean'] 
#' summary(m_aft, par = "log-posterior")[, 'mean']
#' 
#' #----- Frailty model, i.e. site-specific intercepts
#' 
#' m_frail <- stan_surv(
#'   formula = Surv(eventtime, status) ~ trt + (1 | site), 
#'   data    = frail[1:40,], 
#'   basehaz = "exp", 
#'   chains  = 1,
#'   refresh = 0,
#'   iter    = 600,
#'   seed    = 123)
#' print(m_frail)   # shows SD for frailty
#' VarCorr(m_frail) # extract SD explicitly
#'
#' }
#' 
stan_surv <- function(formula, 
                      data, 
                      basehaz          = "ms",
                      basehaz_ops,
                      qnodes           = 15, 
                      prior            = normal(), 
                      prior_intercept  = normal(),
                      prior_aux, 
                      prior_smooth     = exponential(autoscale = FALSE), 
                      prior_covariance = decov(),
                      prior_PD         = FALSE,
                      algorithm        = c("sampling", "meanfield", "fullrank"),
                      adapt_delta      = 0.95, ...) {

  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function.")
  
  if (missing(basehaz_ops)) 
    basehaz_ops <- NULL
  if (missing(data) || !inherits(data, "data.frame"))
    stop("'data' must be a data frame.")
  
  dots      <- list(...)
  algorithm <- match.arg(algorithm)
  
  formula <- parse_formula_and_data(formula, data)
  data    <- formula$data; formula[["data"]] <- NULL
  
  #----------------
  # Construct data
  #----------------
 
  #----- model frame stuff
  
  mf_stuff <- make_model_frame(formula$tf_form, data, drop.unused.levels = TRUE)
  
  mf <- mf_stuff$mf # model frame
  mt <- mf_stuff$mt # model terms
  
  #----- dimensions and response vectors
  
  # entry and exit times for each row of data
  t_beg <- make_t(mf, type = "beg") # entry time
  t_end <- make_t(mf, type = "end") # exit  time
  t_upp <- make_t(mf, type = "upp") # upper time for interval censoring
  
  # ensure no event or censoring times are zero (leads to degenerate
  # estimate for log hazard for most baseline hazards, due to log(0))
  check1 <- any(t_end <= 0, na.rm = TRUE)
  check2 <- any(t_upp <= 0, na.rm = TRUE)
  if (check1 || check2)
    stop2("All event and censoring times must be greater than 0.")

  # event indicator for each row of data
  status <- make_d(mf)
  
  if (any(is.na(status)))
    stop2("Invalid status indicator in Surv object.")
  
  if (any(status < 0 || status > 3))
    stop2("Invalid status indicator in Surv object.")

  # delayed entry indicator for each row of data
  delayed  <- as.logical(!t_beg == 0)
  
  # time variables for stan
  t_event <- aa(t_end[status == 1]) # exact event time
  t_lcens <- aa(t_end[status == 2]) # left  censoring time
  t_rcens <- aa(t_end[status == 0]) # right censoring time
  t_icenl <- aa(t_end[status == 3]) # lower limit of interval censoring time
  t_icenu <- aa(t_upp[status == 3]) # upper limit of interval censoring time
  t_delay <- aa(t_beg[delayed])     # delayed entry time

  # calculate log crude event rate
  t_tmp <- sum(rowMeans(cbind(t_end, t_upp), na.rm = TRUE) - t_beg)
  d_tmp <- sum(!status == 0)
  log_crude_event_rate <- log(d_tmp / t_tmp)
  if (is.infinite(log_crude_event_rate))
    log_crude_event_rate <- 0 # avoids error when there are zero events
  
  # dimensions
  nevent <- sum(status == 1)
  nrcens <- sum(status == 0)
  nlcens <- sum(status == 2)
  nicens <- sum(status == 3)
  ndelay <- sum(delayed)
  
  #----- baseline hazard

  ok_basehaz <- c("exp",
                  "exp-aft",
                  "weibull",
                  "weibull-aft",
                  "gompertz", 
                  "ms", 
                  "bs")
  basehaz <- handle_basehaz_surv(basehaz        = basehaz, 
                                 basehaz_ops    = basehaz_ops,
                                 ok_basehaz     = ok_basehaz,
                                 times          = t_end, 
                                 status         = status,
                                 min_t          = min(t_beg),
                                 max_t          = max(c(t_end,t_upp), na.rm = TRUE))
  nvars <- basehaz$nvars # number of basehaz aux parameters
  
  # flag if intercept is required for baseline hazard
  has_intercept <- ai(has_intercept(basehaz))

  # flag if AFT specification
  is_aft <- get_basehaz_name(basehaz) %in% c("exp-aft", "weibull-aft")
  
  #----- define dimensions and times for quadrature

  # flag if formula uses time-varying effects
  has_tve <- !is.null(formula$td_form)

  # flag if closed form available for cumulative baseline hazard
  has_closed_form <- check_for_closed_form(basehaz)

  # flag for quadrature
  has_quadrature <- has_tve || !has_closed_form
  
  if (has_quadrature) { # model uses quadrature
    
    # standardised nodes and weights for quadrature
    qq <- get_quadpoints(nodes = qnodes)
    qp <- qq$points
    qw <- qq$weights
    
    # quadrature points, evaluated for each row of data
    qpts_event <- uapply(qp, unstandardise_qpts, 0, t_event)
    qpts_lcens <- uapply(qp, unstandardise_qpts, 0, t_lcens)
    qpts_rcens <- uapply(qp, unstandardise_qpts, 0, t_rcens)
    qpts_icenl <- uapply(qp, unstandardise_qpts, 0, t_icenl)
    qpts_icenu <- uapply(qp, unstandardise_qpts, 0, t_icenu)
    qpts_delay <- uapply(qp, unstandardise_qpts, 0, t_delay)

    # quadrature weights, evaluated for each row of data
    qwts_event <- uapply(qw, unstandardise_qwts, 0, t_event)
    qwts_lcens <- uapply(qw, unstandardise_qwts, 0, t_lcens)
    qwts_rcens <- uapply(qw, unstandardise_qwts, 0, t_rcens)
    qwts_icenl <- uapply(qw, unstandardise_qwts, 0, t_icenl)
    qwts_icenu <- uapply(qw, unstandardise_qwts, 0, t_icenu)
    qwts_delay <- uapply(qw, unstandardise_qwts, 0, t_delay)
    
    # times at events and all quadrature points
    cpts_list <- list(t_event,
                      qpts_event,
                      qpts_lcens,
                      qpts_rcens,
                      qpts_icenl,
                      qpts_icenu,
                      qpts_delay)
    idx_cpts <- get_idx_array(sapply(cpts_list, length))
    cpts     <- unlist(cpts_list) # as vector 
    
    # number of quadrature points
    qevent <- length(qwts_event)
    qlcens <- length(qwts_lcens)
    qrcens <- length(qwts_rcens)
    qicens <- length(qwts_icenl)
    qdelay <- length(qwts_delay)

  } else {
    
    # times at all different event types
    cpts_list <- list(t_event,
                      t_lcens,
                      t_rcens,
                      t_icenl,
                      t_icenu,
                      t_delay)
    idx_cpts <- get_idx_array(sapply(cpts_list, length))    
    cpts     <- unlist(cpts_list) # as vector
    
    # dud entries for stan
    qpts_event <- rep(0,0) 
    qpts_lcens <- rep(0,0)
    qpts_rcens <- rep(0,0)
    qpts_icenl <- rep(0,0)
    qpts_icenu <- rep(0,0)
    qpts_delay <- rep(0,0)
    
    if (!qnodes == 15) # warn user if qnodes is not equal to the default
      warning2("There is no quadrature required so 'qnodes' is being ignored.")
    
  }

  #----- basis terms for baseline hazard

  if (!has_quadrature) {
    
    basis_event  <- make_basis(t_event, basehaz)
    ibasis_event <- make_basis(t_event, basehaz, integrate = TRUE)
    ibasis_lcens <- make_basis(t_lcens, basehaz, integrate = TRUE)
    ibasis_rcens <- make_basis(t_rcens, basehaz, integrate = TRUE)
    ibasis_icenl <- make_basis(t_icenl, basehaz, integrate = TRUE)
    ibasis_icenu <- make_basis(t_icenu, basehaz, integrate = TRUE)
    ibasis_delay <- make_basis(t_delay, basehaz, integrate = TRUE)
    
  } else {
    
    basis_epts_event <- make_basis(t_event,    basehaz)
    basis_qpts_event <- make_basis(qpts_event, basehaz)
    basis_qpts_lcens <- make_basis(qpts_lcens, basehaz)
    basis_qpts_rcens <- make_basis(qpts_rcens, basehaz)
    basis_qpts_icenl <- make_basis(qpts_icenl, basehaz)
    basis_qpts_icenu <- make_basis(qpts_icenu, basehaz)
    basis_qpts_delay <- make_basis(qpts_delay, basehaz)
    
  }
    
  #----- model frames for generating predictor matrices

  mf_event <- keep_rows(mf, status == 1)
  mf_lcens <- keep_rows(mf, status == 2)
  mf_rcens <- keep_rows(mf, status == 0)
  mf_icens <- keep_rows(mf, status == 3)
  mf_delay <- keep_rows(mf, delayed)  
   
  if (!has_quadrature) {
    
    # combined model frame, without quadrature
    mf_cpts <- rbind(mf_event,
                     mf_lcens,
                     mf_rcens,
                     mf_icens,
                     mf_icens,
                     mf_delay)
    
  } else {
    
    # combined model frame, with quadrature
    mf_cpts <- rbind(mf_event,
                     rep_rows(mf_event, times = qnodes),
                     rep_rows(mf_lcens, times = qnodes),
                     rep_rows(mf_rcens, times = qnodes),
                     rep_rows(mf_icens, times = qnodes),
                     rep_rows(mf_icens, times = qnodes),
                     rep_rows(mf_delay, times = qnodes))

  }

  if (has_tve) {
    
    # generate a model frame with time transformations for tve effects
    mf_tve <- make_model_frame(formula$tt_frame, data.frame(times__ = cpts))$mf
    
    # NB next line avoids dropping terms attribute from 'mf_cpts'
    mf_cpts[, colnames(mf_tve)] <- mf_tve
    
  }  

  #----- time-fixed predictor matrices
  
  ff        <- formula$fe_form
  x         <- make_x(ff, mf     )$x
  x_cpts    <- make_x(ff, mf_cpts)$x
  x_centred <- sweep(x_cpts, 2, colMeans(x), FUN = "-")
  K         <- ncol(x_cpts)
  
  if (!has_quadrature) {
    
    # time-fixed predictor matrices, without quadrature
    # NB skip index 5 on purpose, since time fixed predictor matrix is 
    # identical for lower and upper limits of interval censoring time
    x_event <- x_centred[idx_cpts[1,1]:idx_cpts[1,2], , drop = FALSE]
    x_lcens <- x_centred[idx_cpts[2,1]:idx_cpts[2,2], , drop = FALSE]
    x_rcens <- x_centred[idx_cpts[3,1]:idx_cpts[3,2], , drop = FALSE]
    x_icens <- x_centred[idx_cpts[4,1]:idx_cpts[4,2], , drop = FALSE]
    x_delay <- x_centred[idx_cpts[6,1]:idx_cpts[6,2], , drop = FALSE]
    
  } else {
    
    # time-fixed predictor matrices, with quadrature
    # NB skip index 6 on purpose, since time fixed predictor matrix is 
    # identical for lower and upper limits of interval censoring time
    x_epts_event <- x_centred[idx_cpts[1,1]:idx_cpts[1,2], , drop = FALSE]
    x_qpts_event <- x_centred[idx_cpts[2,1]:idx_cpts[2,2], , drop = FALSE]
    x_qpts_lcens <- x_centred[idx_cpts[3,1]:idx_cpts[3,2], , drop = FALSE]
    x_qpts_rcens <- x_centred[idx_cpts[4,1]:idx_cpts[4,2], , drop = FALSE]
    x_qpts_icens <- x_centred[idx_cpts[5,1]:idx_cpts[5,2], , drop = FALSE]
    x_qpts_delay <- x_centred[idx_cpts[7,1]:idx_cpts[7,2], , drop = FALSE]
 
  }
  
  #----- time-varying predictor matrices
  
  if (has_tve) {
    
    # time-varying predictor matrix
    s_cpts          <- make_s(formula, mf_cpts, xlevs = xlevs)
    smooth_map      <- get_smooth_name(s_cpts, type = "smooth_map")
    smooth_idx      <- get_idx_array(table(smooth_map))
    S <- ncol(s_cpts) # number of tve coefficients
    
    # store some additional information in model formula
    # stating how many columns in the predictor matrix
    # each tve() term in the model formula corresponds to
    formula$tt_ncol <- attr(s_cpts, "tt_ncol")
    formula$tt_map  <- attr(s_cpts, "tt_map")
    
  } else {
    
    # dud entries if no tve() terms in model formula
    s_cpts          <- matrix(0,length(cpts),0)
    smooth_idx      <- matrix(0,0,2)
    smooth_map      <- integer(0)
    S               <- 0L
    
    formula$tt_ncol <- integer(0)
    formula$tt_map  <- integer(0)
    
  }
  
  if (has_quadrature) {
    
    # time-varying predictor matrices, with quadrature
    s_epts_event <- s_cpts[idx_cpts[1,1]:idx_cpts[1,2], , drop = FALSE]
    s_qpts_event <- s_cpts[idx_cpts[2,1]:idx_cpts[2,2], , drop = FALSE]
    s_qpts_lcens <- s_cpts[idx_cpts[3,1]:idx_cpts[3,2], , drop = FALSE]
    s_qpts_rcens <- s_cpts[idx_cpts[4,1]:idx_cpts[4,2], , drop = FALSE]
    s_qpts_icenl <- s_cpts[idx_cpts[5,1]:idx_cpts[5,2], , drop = FALSE]
    s_qpts_icenu <- s_cpts[idx_cpts[6,1]:idx_cpts[6,2], , drop = FALSE]
    s_qpts_delay <- s_cpts[idx_cpts[7,1]:idx_cpts[7,2], , drop = FALSE]
  
  }

  #----- random effects predictor matrices

  has_bars <- as.logical(length(formula$bars))
  
  # use 'stan_glmer' approach
  if (has_bars) {
    
    group_unpadded <- lme4::mkReTrms(formula$bars, mf_cpts)  
    group <- pad_reTrms(Ztlist = group_unpadded$Ztlist,
                        cnms   = group_unpadded$cnms,
                        flist  = group_unpadded$flist)
    z_cpts <- group$Z
    
  } else {
    
    group  <- NULL
    z_cpts <- matrix(0,length(cpts),0)
    
  } 
  
  if (!has_quadrature) {
    
    # random effects predictor matrices, without quadrature
    # NB skip index 5 on purpose, since time fixed predictor matrix is 
    # identical for lower and upper limits of interval censoring time
    z_event <- z_cpts[idx_cpts[1,1]:idx_cpts[1,2], , drop = FALSE]
    z_lcens <- z_cpts[idx_cpts[2,1]:idx_cpts[2,2], , drop = FALSE]
    z_rcens <- z_cpts[idx_cpts[3,1]:idx_cpts[3,2], , drop = FALSE]
    z_icens <- z_cpts[idx_cpts[4,1]:idx_cpts[4,2], , drop = FALSE]
    z_delay <- z_cpts[idx_cpts[6,1]:idx_cpts[6,2], , drop = FALSE]
    
    parts_event <- extract_sparse_parts(z_event)
    parts_lcens <- extract_sparse_parts(z_lcens)
    parts_rcens <- extract_sparse_parts(z_rcens)
    parts_icens <- extract_sparse_parts(z_icens)
    parts_delay <- extract_sparse_parts(z_delay)
    
  } else {
    
    # random effects predictor matrices, with quadrature
    # NB skip index 6 on purpose, since time fixed predictor matrix is 
    # identical for lower and upper limits of interval censoring time
    z_epts_event <- z_cpts[idx_cpts[1,1]:idx_cpts[1,2], , drop = FALSE]
    z_qpts_event <- z_cpts[idx_cpts[2,1]:idx_cpts[2,2], , drop = FALSE]
    z_qpts_lcens <- z_cpts[idx_cpts[3,1]:idx_cpts[3,2], , drop = FALSE]
    z_qpts_rcens <- z_cpts[idx_cpts[4,1]:idx_cpts[4,2], , drop = FALSE]
    z_qpts_icens <- z_cpts[idx_cpts[5,1]:idx_cpts[5,2], , drop = FALSE]
    z_qpts_delay <- z_cpts[idx_cpts[7,1]:idx_cpts[7,2], , drop = FALSE]
    
    parts_epts_event <- extract_sparse_parts(z_epts_event)
    parts_qpts_event <- extract_sparse_parts(z_qpts_event)
    parts_qpts_lcens <- extract_sparse_parts(z_qpts_lcens)
    parts_qpts_rcens <- extract_sparse_parts(z_qpts_rcens)
    parts_qpts_icens <- extract_sparse_parts(z_qpts_icens)
    parts_qpts_delay <- extract_sparse_parts(z_qpts_delay)
    
  }

  #----- stan data
  
  standata <- nlist(
    K, S, 
    nvars,
    x_bar = aa(colMeans(x)),
    has_intercept, 
    has_quadrature,
    smooth_map,
    smooth_idx,
    type = basehaz$type,
    log_crude_event_rate = 
      ifelse(is_aft, -log_crude_event_rate, log_crude_event_rate),

    nevent       = if (has_quadrature) 0L else nevent,
    nlcens       = if (has_quadrature) 0L else nlcens,
    nrcens       = if (has_quadrature) 0L else nrcens,
    nicens       = if (has_quadrature) 0L else nicens,
    ndelay       = if (has_quadrature) 0L else ndelay,
    
    t_event      = if (has_quadrature) rep(0,0) else t_event,
    t_lcens      = if (has_quadrature) rep(0,0) else t_lcens,
    t_rcens      = if (has_quadrature) rep(0,0) else t_rcens,
    t_icenl      = if (has_quadrature) rep(0,0) else t_icenl,
    t_icenu      = if (has_quadrature) rep(0,0) else t_icenu,
    t_delay      = if (has_quadrature) rep(0,0) else t_delay,
    
    x_event      = if (has_quadrature) matrix(0,0,K) else x_event,
    x_lcens      = if (has_quadrature) matrix(0,0,K) else x_lcens,
    x_rcens      = if (has_quadrature) matrix(0,0,K) else x_rcens,
    x_icens      = if (has_quadrature) matrix(0,0,K) else x_icens,
    x_delay      = if (has_quadrature) matrix(0,0,K) else x_delay,

    w_event      = if (has_quadrature || !has_bars || nevent == 0) double(0) else parts_event$w,
    w_lcens      = if (has_quadrature || !has_bars || nlcens == 0) double(0) else parts_lcens$w,
    w_rcens      = if (has_quadrature || !has_bars || nrcens == 0) double(0) else parts_rcens$w,
    w_icens      = if (has_quadrature || !has_bars || nicens == 0) double(0) else parts_icens$w,
    w_delay      = if (has_quadrature || !has_bars || ndelay == 0) double(0) else parts_delay$w,
 
    v_event      = if (has_quadrature || !has_bars || nevent == 0) integer(0) else parts_event$v - 1L,
    v_lcens      = if (has_quadrature || !has_bars || nlcens == 0) integer(0) else parts_lcens$v - 1L,
    v_rcens      = if (has_quadrature || !has_bars || nrcens == 0) integer(0) else parts_rcens$v - 1L,
    v_icens      = if (has_quadrature || !has_bars || nicens == 0) integer(0) else parts_icens$v - 1L,
    v_delay      = if (has_quadrature || !has_bars || ndelay == 0) integer(0) else parts_delay$v - 1L,    
    
    u_event      = if (has_quadrature || !has_bars || nevent == 0) integer(0) else parts_event$u - 1L,
    u_lcens      = if (has_quadrature || !has_bars || nlcens == 0) integer(0) else parts_lcens$u - 1L,
    u_rcens      = if (has_quadrature || !has_bars || nrcens == 0) integer(0) else parts_rcens$u - 1L,
    u_icens      = if (has_quadrature || !has_bars || nicens == 0) integer(0) else parts_icens$u - 1L,
    u_delay      = if (has_quadrature || !has_bars || ndelay == 0) integer(0) else parts_delay$u - 1L,    
  
    nnz_event    = if (has_quadrature || !has_bars || nevent == 0) 0L else length(parts_event$w),
    nnz_lcens    = if (has_quadrature || !has_bars || nlcens == 0) 0L else length(parts_lcens$w),
    nnz_rcens    = if (has_quadrature || !has_bars || nrcens == 0) 0L else length(parts_rcens$w),
    nnz_icens    = if (has_quadrature || !has_bars || nicens == 0) 0L else length(parts_icens$w),
    nnz_delay    = if (has_quadrature || !has_bars || ndelay == 0) 0L else length(parts_delay$w),    
    
    basis_event  = if (has_quadrature) matrix(0,0,nvars) else basis_event,
    ibasis_event = if (has_quadrature) matrix(0,0,nvars) else ibasis_event,
    ibasis_lcens = if (has_quadrature) matrix(0,0,nvars) else ibasis_lcens,
    ibasis_rcens = if (has_quadrature) matrix(0,0,nvars) else ibasis_rcens,
    ibasis_icenl = if (has_quadrature) matrix(0,0,nvars) else ibasis_icenl,
    ibasis_icenu = if (has_quadrature) matrix(0,0,nvars) else ibasis_icenu,
    ibasis_delay = if (has_quadrature) matrix(0,0,nvars) else ibasis_delay,

    qnodes       = if (!has_quadrature) 0L else qnodes,
    
    Nevent       = if (!has_quadrature) 0L else nevent,
    Nlcens       = if (!has_quadrature) 0L else nlcens, 
    Nrcens       = if (!has_quadrature) 0L else nrcens, 
    Nicens       = if (!has_quadrature) 0L else nicens,
    Ndelay       = if (!has_quadrature) 0L else ndelay,
    
    qevent       = if (!has_quadrature) 0L else qevent,
    qlcens       = if (!has_quadrature) 0L else qlcens,
    qrcens       = if (!has_quadrature) 0L else qrcens,
    qicens       = if (!has_quadrature) 0L else qicens,
    qdelay       = if (!has_quadrature) 0L else qdelay,
    
    epts_event   = if (!has_quadrature) rep(0,0) else t_event,
    qpts_event   = if (!has_quadrature) rep(0,0) else qpts_event,
    qpts_lcens   = if (!has_quadrature) rep(0,0) else qpts_lcens,
    qpts_rcens   = if (!has_quadrature) rep(0,0) else qpts_rcens,
    qpts_icenl   = if (!has_quadrature) rep(0,0) else qpts_icenl,
    qpts_icenu   = if (!has_quadrature) rep(0,0) else qpts_icenu,
    qpts_delay   = if (!has_quadrature) rep(0,0) else qpts_delay,
    
    qwts_event   = if (!has_quadrature) rep(0,0) else qwts_event,
    qwts_lcens   = if (!has_quadrature) rep(0,0) else qwts_lcens,
    qwts_rcens   = if (!has_quadrature) rep(0,0) else qwts_rcens,
    qwts_icenl   = if (!has_quadrature) rep(0,0) else qwts_icenl,
    qwts_icenu   = if (!has_quadrature) rep(0,0) else qwts_icenu,
    qwts_delay   = if (!has_quadrature) rep(0,0) else qwts_delay,
    
    x_epts_event = if (!has_quadrature) matrix(0,0,K) else x_epts_event,
    x_qpts_event = if (!has_quadrature) matrix(0,0,K) else x_qpts_event,
    x_qpts_lcens = if (!has_quadrature) matrix(0,0,K) else x_qpts_lcens,
    x_qpts_rcens = if (!has_quadrature) matrix(0,0,K) else x_qpts_rcens,
    x_qpts_icens = if (!has_quadrature) matrix(0,0,K) else x_qpts_icens,
    x_qpts_delay = if (!has_quadrature) matrix(0,0,K) else x_qpts_delay,
    
    s_epts_event = if (!has_quadrature) matrix(0,0,S) else s_epts_event,
    s_qpts_event = if (!has_quadrature) matrix(0,0,S) else s_qpts_event,
    s_qpts_lcens = if (!has_quadrature) matrix(0,0,S) else s_qpts_lcens,
    s_qpts_rcens = if (!has_quadrature) matrix(0,0,S) else s_qpts_rcens,
    s_qpts_icenl = if (!has_quadrature) matrix(0,0,S) else s_qpts_icenl,
    s_qpts_icenu = if (!has_quadrature) matrix(0,0,S) else s_qpts_icenu,
    s_qpts_delay = if (!has_quadrature) matrix(0,0,S) else s_qpts_delay,
    
    w_epts_event = if (!has_quadrature || !has_bars || qevent == 0) double(0) else parts_epts_event$w,
    w_qpts_event = if (!has_quadrature || !has_bars || qevent == 0) double(0) else parts_qpts_event$w,
    w_qpts_lcens = if (!has_quadrature || !has_bars || qlcens == 0) double(0) else parts_qpts_lcens$w,
    w_qpts_rcens = if (!has_quadrature || !has_bars || qrcens == 0) double(0) else parts_qpts_rcens$w,
    w_qpts_icens = if (!has_quadrature || !has_bars || qicens == 0) double(0) else parts_qpts_icens$w,
    w_qpts_delay = if (!has_quadrature || !has_bars || qdelay == 0) double(0) else parts_qpts_delay$w,
    
    v_epts_event = if (!has_quadrature || !has_bars || qevent == 0) integer(0) else parts_epts_event$v - 1L,
    v_qpts_event = if (!has_quadrature || !has_bars || qevent == 0) integer(0) else parts_qpts_event$v - 1L,
    v_qpts_lcens = if (!has_quadrature || !has_bars || qlcens == 0) integer(0) else parts_qpts_lcens$v - 1L,
    v_qpts_rcens = if (!has_quadrature || !has_bars || qrcens == 0) integer(0) else parts_qpts_rcens$v - 1L,
    v_qpts_icens = if (!has_quadrature || !has_bars || qicens == 0) integer(0) else parts_qpts_icens$v - 1L,
    v_qpts_delay = if (!has_quadrature || !has_bars || qdelay == 0) integer(0) else parts_qpts_delay$v - 1L,    
    
    u_epts_event = if (!has_quadrature || !has_bars || qevent == 0) integer(0) else parts_epts_event$u - 1L,
    u_qpts_event = if (!has_quadrature || !has_bars || qevent == 0) integer(0) else parts_qpts_event$u - 1L,
    u_qpts_lcens = if (!has_quadrature || !has_bars || qlcens == 0) integer(0) else parts_qpts_lcens$u - 1L,
    u_qpts_rcens = if (!has_quadrature || !has_bars || qrcens == 0) integer(0) else parts_qpts_rcens$u - 1L,
    u_qpts_icens = if (!has_quadrature || !has_bars || qicens == 0) integer(0) else parts_qpts_icens$u - 1L,
    u_qpts_delay = if (!has_quadrature || !has_bars || qdelay == 0) integer(0) else parts_qpts_delay$u - 1L,    
    
    nnz_epts_event = if (!has_quadrature || !has_bars || qevent == 0) 0L else length(parts_epts_event$w),
    nnz_qpts_event = if (!has_quadrature || !has_bars || qevent == 0) 0L else length(parts_qpts_event$w),
    nnz_qpts_lcens = if (!has_quadrature || !has_bars || qlcens == 0) 0L else length(parts_qpts_lcens$w),
    nnz_qpts_rcens = if (!has_quadrature || !has_bars || qrcens == 0) 0L else length(parts_qpts_rcens$w),
    nnz_qpts_icens = if (!has_quadrature || !has_bars || qicens == 0) 0L else length(parts_qpts_icens$w),
    nnz_qpts_delay = if (!has_quadrature || !has_bars || qdelay == 0) 0L else length(parts_qpts_delay$w),
    
    basis_epts_event = if (!has_quadrature) matrix(0,0,nvars) else basis_epts_event,
    basis_qpts_event = if (!has_quadrature) matrix(0,0,nvars) else basis_qpts_event,
    basis_qpts_lcens = if (!has_quadrature) matrix(0,0,nvars) else basis_qpts_lcens,
    basis_qpts_rcens = if (!has_quadrature) matrix(0,0,nvars) else basis_qpts_rcens,
    basis_qpts_icenl = if (!has_quadrature) matrix(0,0,nvars) else basis_qpts_icenl,
    basis_qpts_icenu = if (!has_quadrature) matrix(0,0,nvars) else basis_qpts_icenu,
    basis_qpts_delay = if (!has_quadrature) matrix(0,0,nvars) else basis_qpts_delay
  )
 
  #----- random-effects structure
  
  if (has_bars) {
    
    fl <- group$flist
    p  <- sapply(group$cnms, FUN = length)
    l  <- sapply(attr(fl, "assign"), function(i) nlevels(fl[[i]]))
    t  <- length(l)
    standata$p <- as.array(p)   # num ranefs for each grouping factor
    standata$l <- as.array(l)   # num levels for each grouping factor
    standata$t <- t             # num of grouping factors
    standata$q <- ncol(group$Z) # p * l
    standata$special_case <- all(sapply(group$cnms, intercept_only))
    
  } else { # no random effects structure
    
    standata$p <- integer(0)
    standata$l <- integer(0)
    standata$t <- 0L
    standata$q <- 0L
    standata$special_case <- 0L
    
  }

  #----- priors and hyperparameters

  # valid priors
  ok_dists <- nlist("normal", 
                    student_t = "t", 
                    "cauchy", 
                    "hs", 
                    "hs_plus", 
                    "laplace", 
                    "lasso") # disallow product normal
  ok_intercept_dists  <- ok_dists[1:3]
  ok_aux_dists        <- get_ok_priors_for_aux(basehaz)
  ok_smooth_dists     <- c(ok_dists[1:3], "exponential")
  ok_covariance_dists <- c("decov")
  
  if (missing(prior_aux))
    prior_aux <- get_default_prior_for_aux(basehaz)
  
  # priors
  user_prior_stuff <- prior_stuff <-
    handle_glm_prior(prior, 
                     nvars = K,
                     default_scale = 2.5,
                     link = NULL,
                     ok_dists = ok_dists)
  
  user_prior_intercept_stuff <- prior_intercept_stuff <-
    handle_glm_prior(prior_intercept, 
                     nvars = 1,
                     default_scale = 20,
                     link = NULL,
                     ok_dists = ok_intercept_dists)
 
  user_prior_aux_stuff <- prior_aux_stuff <-
    handle_glm_prior(prior_aux, 
                     nvars = basehaz$nvars,
                     default_scale = get_default_aux_scale(basehaz),
                     link = NULL,
                     ok_dists = ok_aux_dists)

  user_prior_smooth_stuff <- prior_smooth_stuff <-
    handle_glm_prior(prior_smooth, 
                     nvars = if (S) max(smooth_map) else 0,
                     default_scale = 1,
                     link = NULL,
                     ok_dists = ok_smooth_dists)
    
  # stop null priors when prior_PD is true
  if (prior_PD) {
    if (is.null(prior))
      stop("'prior' cannot be NULL if 'prior_PD' is TRUE.")
    if (is.null(prior_intercept) && has_intercept)
      stop("'prior_intercept' cannot be NULL if 'prior_PD' is TRUE.")
    if (is.null(prior_aux))
      stop("'prior_aux' cannot be NULL if 'prior_PD' is TRUE.")    
    if (is.null(prior_smooth) && (S > 0))
      stop("'prior_smooth' cannot be NULL if 'prior_PD' is TRUE.")    
  }
  
  # handle prior for random effects structure
  if (has_bars) {
    
    user_prior_b_stuff <- prior_b_stuff <- 
      handle_cov_prior(prior_covariance, 
                       cnms = group$cnms, 
                       ok_dists = ok_covariance_dists)
    
    if (is.null(prior_covariance))
      stop("'prior_covariance' cannot be NULL.")
  
  } else {
    user_prior_b_stuff <- NULL
    prior_b_stuff      <- NULL
    prior_covariance   <- NULL
  }
  
  # autoscaling of priors
  prior_stuff           <- autoscale_prior(prior_stuff, predictors = x)
  prior_intercept_stuff <- autoscale_prior(prior_intercept_stuff)
  prior_aux_stuff       <- autoscale_prior(prior_aux_stuff)
  prior_smooth_stuff    <- autoscale_prior(prior_smooth_stuff)
  
  # priors
  standata$prior_dist              <- prior_stuff$prior_dist
  standata$prior_dist_for_intercept<- prior_intercept_stuff$prior_dist
  standata$prior_dist_for_aux      <- prior_aux_stuff$prior_dist
  standata$prior_dist_for_smooth   <- prior_smooth_stuff$prior_dist
  standata$prior_dist_for_cov      <- prior_b_stuff$prior_dist
  
  # hyperparameters
  standata$prior_mean               <- prior_stuff$prior_mean
  standata$prior_scale              <- prior_stuff$prior_scale
  standata$prior_df                 <- prior_stuff$prior_df
  standata$prior_mean_for_intercept <- c(prior_intercept_stuff$prior_mean)
  standata$prior_scale_for_intercept<- c(prior_intercept_stuff$prior_scale)
  standata$prior_df_for_intercept   <- c(prior_intercept_stuff$prior_df)
  standata$prior_scale_for_aux      <- prior_aux_stuff$prior_scale
  standata$prior_df_for_aux         <- prior_aux_stuff$prior_df
  standata$prior_conc_for_aux       <- prior_aux_stuff$prior_concentration
  standata$prior_mean_for_smooth    <- prior_smooth_stuff$prior_mean
  standata$prior_scale_for_smooth   <- prior_smooth_stuff$prior_scale
  standata$prior_df_for_smooth      <- prior_smooth_stuff$prior_df
  standata$global_prior_scale       <- prior_stuff$global_prior_scale
  standata$global_prior_df          <- prior_stuff$global_prior_df
  standata$slab_df                  <- prior_stuff$slab_df
  standata$slab_scale               <- prior_stuff$slab_scale
  
  # hyperparameters for covariance
  if (has_bars) {
    standata$b_prior_shape          <- prior_b_stuff$prior_shape
    standata$b_prior_scale          <- prior_b_stuff$prior_scale
    standata$concentration          <- prior_b_stuff$prior_concentration
    standata$regularization         <- prior_b_stuff$prior_regularization
    standata$len_concentration      <- length(standata$concentration)
    standata$len_regularization     <- length(standata$regularization)
    standata$len_theta_L            <- sum(choose(standata$p, 2), standata$p)  
  } else { # no random effects structure
    standata$b_prior_shape          <- rep(0, 0)
    standata$b_prior_scale          <- rep(0, 0)
    standata$concentration          <- rep(0, 0)
    standata$regularization         <- rep(0, 0)
    standata$len_concentration      <- 0L
    standata$len_regularization     <- 0L    
    standata$len_theta_L            <- 0L
  }
  
  # any additional flags
  standata$prior_PD <- ai(prior_PD)

  #---------------
  # Prior summary
  #---------------
  
  prior_info <- summarize_jm_prior(
    user_priorEvent           = user_prior_stuff,
    user_priorEvent_intercept = user_prior_intercept_stuff,
    user_priorEvent_aux       = user_prior_aux_stuff,
    adjusted_priorEvent_scale           = prior_stuff$prior_scale,
    adjusted_priorEvent_intercept_scale = prior_intercept_stuff$prior_scale,
    adjusted_priorEvent_aux_scale       = prior_aux_stuff$prior_scale,
    e_has_intercept  = has_intercept,
    e_has_predictors = K > 0,
    basehaz = basehaz,
    user_prior_covariance = prior_covariance,
    b_user_prior_stuff = user_prior_b_stuff,
    b_prior_stuff = prior_b_stuff
  )
  
  #-----------
  # Fit model
  #-----------

  # obtain stan model code
  stanfit  <- stanmodels$surv
  
  # specify parameters for stan to monitor
  stanpars <- c(if (standata$has_intercept) "alpha",
                if (standata$K)             "beta",
                if (standata$S)             "beta_tve",
                if (standata$S)             "smooth_sd",
                if (standata$nvars)         "aux",
                if (standata$t)             "b",
                if (standata$t)             "theta_L")
  
  # fit model using stan
  if (algorithm == "sampling") { # mcmc
    args <- set_sampling_args(
      object = stanfit, 
      data   = standata, 
      pars   = stanpars, 
      prior  = prior, 
      user_dots = list(...), 
      user_adapt_delta = adapt_delta, 
      show_messages = FALSE)
    stanfit <- do.call(rstan::sampling, args)
  } else { # meanfield or fullrank vb
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      algorithm
    )
    args[names(dots)] <- dots
    stanfit <- do.call(rstan::vb, args)
  }
  check_stanfit(stanfit)
  
  # replace 'theta_L' with the variance-covariance matrix
  if (has_bars)
    stanfit <- evaluate_Sigma(stanfit, group$cnms)
  
  # define new parameter names
  nms_beta   <- colnames(x_cpts) # may be NULL
  nms_tve    <- get_smooth_name(s_cpts, type = "smooth_coefs") # may be NULL
  nms_smooth <- get_smooth_name(s_cpts, type = "smooth_sd")    # may be NULL
  nms_int    <- get_int_name_basehaz(basehaz)
  nms_aux    <- get_aux_name_basehaz(basehaz)
  nms_b      <- get_b_names(group)                             # may be NULL
  nms_vc     <- get_varcov_names(group)                        # may be NULL
  nms_all    <- c(nms_int,
                  nms_beta,
                  nms_tve,
                  nms_smooth,
                  nms_aux,
                  nms_b,
                  nms_vc,
                  "log-posterior")

  # substitute new parameter names into 'stanfit' object
  stanfit <- replace_stanfit_nms(stanfit, nms_all)
  
  # return an object of class 'stansurv'
  fit <- nlist(stanfit, 
               formula,
               has_tve,
               has_quadrature,
               has_bars,
               data,
               model_frame      = mf,
               terms            = mt,
               xlevels          = .getXlevels(mt, mf),
               x,
               x_cpts,
               s_cpts           = if (has_tve)  s_cpts else NULL,
               z_cpts           = if (has_bars) z_cpts else NULL,
               cnms             = if (has_bars) group_unpadded$cnms  else NULL,
               flist            = if (has_bars) group_unpadded$flist else NULL,
               t_beg, 
               t_end,
               status,
               event            = as.logical(status == 1),
               delayed,
               basehaz,
               nobs             = nrow(mf),
               nevents          = nevent,
               nlcens,
               nrcens,
               nicens,
               ncensor          = nlcens + nrcens + nicens,
               ndelayed         = ndelay,
               prior_info,
               qnodes           = if (has_quadrature) qnodes else NULL,
               algorithm,
               stan_function    = "stan_surv",
               rstanarm_version = utils::packageVersion("rstanarm"),
               call             = match.call(expand.dots = TRUE))
  stansurv(fit)
}


#' Time-varying effects in Bayesian survival models
#' 
#' This is a special function that can be used in the formula of a Bayesian 
#' survival model estimated using \code{\link{stan_surv}}. It specifies that a 
#' time-varying coefficient should be estimated for the covariate \code{x}. 
#' The time-varying coefficient is currently modelled using B-splines (with
#' piecewise constant included as a special case). Note that the \code{tve} 
#' function only has meaning when evaluated within the formula of a 
#' \code{\link{stan_surv}} call and does not have meaning outside of that 
#' context. The exported function documented here just returns \code{x}. 
#' However when called internally the \code{tve} function returns several 
#' other pieces of useful information used in the model fitting.
#' 
#' @export
#' 
#' @param x The covariate for which a time-varying coefficient should be 
#'   estimated.
#' @param type The type of function used to model the time-varying coefficient.
#'   Currently only \code{type = "bs"} is allowed. This corresponds to a
#'   B-spline function. Note that \emph{cubic} B-splines are used by default
#'   but this can be changed by the user via the \code{degree} argument 
#'   described below. Of particular note is that \code{degree = 0} is 
#'   is treated as a special case corresponding to a piecewise constant basis.
#' @param df A positive integer specifying the degrees of freedom 
#'   for the B-spline function. Two boundary knots and \code{df - degree} 
#'   internal knots are used to generate the B-spline function.
#'   The internal knots are placed at equally spaced percentiles of the 
#'   distribution of the uncensored event times. The default is to use 
#'   \code{df = 3} unless \code{df} or \code{knots} is explicitly 
#'   specified by the user.
#' @param knots A numeric vector explicitly specifying internal knot  
#'   locations for the B-spline function. Note that \code{knots} cannot be 
#'   specified if \code{df} is specified. Also note that this argument only 
#'   controls the \emph{internal} knot locations. In addition, boundary 
#'   knots are placed at the earliest entry time and latest event or 
#'   censoring time and these cannot be changed by the user.
#' @param degree A positive integer specifying the degree for the B-spline 
#'   function. The order of the B-spline is equal to \code{degree + 1}.
#'   Note that \code{degree = 0} is allowed and is treated as a special
#'   case corresponding to a piecewise constant basis.
#'     
#' @return The exported \code{tve} function documented here just returns 
#'   \code{x}. However, when called internally the \code{tve} function returns 
#'   several other pieces of useful information. For the most part, these are 
#'   added to the formula element of the returned \code{\link{stanreg}} object
#'   (that is \code{object[["formula"]]} where \code{object} is the fitted
#'   model). Information added to the formula element of the \code{stanreg} 
#'   object includes the following:
#'   \itemize{
#'   \item \code{tt_vars}: A list with the names of variables in the model
#'   formula that were wrapped in the \code{tve} function.
#'   \item \code{tt_types}: A list with the \code{type} (e.g. \code{"bs"}) 
#'   of \code{tve} function corresponding to each variable in \code{tt_vars}.
#'   \item \code{tt_degrees}: A list with the \code{degree} for the
#'   B-spline function corresponding to each variable in \code{tt_vars}.
#'   \item \code{tt_calls}: A list with the call required to construct the
#'   transformation of time for each variable in \code{tt_vars}.
#'   \item \code{tt_forms}: Same as \code{tt_calls} but expressed as formulas.
#'   \item \code{tt_frame}: A single formula that can be used to generate a 
#'   model frame that contains the unique set of transformations of time 
#'   (i.e. the basis terms) that are required to build all time-varying 
#'   coefficients in the model. In other words a single formula with the 
#'   unique element(s) contained in \code{tt_forms}.
#'   }
#'
#' @examples 
#' # Exported function just returns the input variable
#' identical(pbcSurv$trt, tve(pbcSurv$trt)) # returns TRUE
#' 
#' # Internally the function returns and stores information 
#' # used to form the time-varying coefficients in the model
#' m1 <- stan_surv(Surv(futimeYears, death) ~ tve(trt) + tve(sex, degree = 0),
#'                 data = pbcSurv, chains = 1, iter = 50)
#' m1$formula[["tt_vars"]]
#' m1$formula[["tt_forms"]]
#'  
tve <- function(x, 
                type   = "bs", 
                df     = NULL, 
                knots  = NULL,
                degree = 3L) {
  
  type <- match.arg(type)
  
  if (!is.null(df) && !is.null(knots))
    stop("Cannot specify both 'df' and 'knots' in the 'tve' function.")
  
  if (degree < 0)
    stop("In 'tve' function, 'degree' must be non-negative.")
  
  if (is.null(df) && is.null(knots))
    df <- 3L
  
  x
}


#---------- internal

# Construct a list with information about the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements: df, knots, degree
# @param ok_basehaz A list of admissible baseline hazards
# @param times A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @param min_t Scalar, the minimum entry time across all individuals
# @param max_t Scalar, the maximum event or censoring time across all individuals
# @return A named list with the following elements:
#   type: integer specifying the type of baseline hazard, 1L = weibull,
#     2L = b-splines, 3L = piecewise.
#   type_name: character string specifying the type of baseline hazard.
#   user_df: integer specifying the input to the df argument
#   df: integer specifying the number of parameters to use for the 
#     baseline hazard.
#   knots: the knot locations for the baseline hazard.
#   bs_basis: The basis terms for the B-splines. This is passed to Stan
#     as the "model matrix" for the baseline hazard. It is also used in
#     post-estimation when evaluating the baseline hazard for posterior
#     predictions since it contains information about the knot locations
#     for the baseline hazard (this is implemented via splines::predict.bs). 
handle_basehaz_surv <- function(basehaz, 
                                basehaz_ops,
                                ok_basehaz,
                                times, 
                                status,
                                min_t, max_t) {
  
  if (!basehaz %in% ok_basehaz)
    stop2("'basehaz' should be one of: ", comma(ok_basehaz))
  
  ok_basehaz_ops <- get_ok_basehaz_ops(basehaz)
  if (!all(names(basehaz_ops) %in% ok_basehaz_ops))
    stop2("'basehaz_ops' can only include: ", comma(ok_basehaz_ops))

  if (basehaz %in% c("ms", "bs", "piecewise")) {
    
    df     <- basehaz_ops$df
    knots  <- basehaz_ops$knots
    degree <- basehaz_ops$degree
    
    if (!is.null(df) && !is.null(knots))
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    
    if (is.null(df))
      df <- switch(basehaz,
                   "ms"        = 6L, # assumes intercept
                   "bs"        = 5L, # assumes no intercept
                   "piecewise" = 5L, # assumes no intercept
                   df) # NB this is ignored if the user specified knots
    
    if (is.null(degree))
      degree <- 3L # cubic splines
    
    tt <- times[status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
      warning2("No observed events found in the data. Censoring times will ",
               "be used to evaluate default knot locations for splines.")
      tt <- times
    }
    
    if (!is.null(knots)) {
      if (any(knots < min_t))
        stop2("'knots' cannot be placed before the earliest entry time.")
      if (any(knots > max_t))
        stop2("'knots' cannot be placed beyond the latest event time.")
    }
    
  }
  
  if (basehaz %in% c("exp", "exp-aft")) {
    
    degree <- NULL # degree for splines
    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 0L   # number of aux parameters, none
    
  } else if (basehaz %in% c("weibull", "weibull-aft")) {
    
    degree <- NULL # degree for splines
    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 1L   # number of aux parameters, Weibull shape
    
  } else if (basehaz == "gompertz") {
    
    degree <- NULL # degree for splines
    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 1L   # number of aux parameters, Gompertz scale
    
  } else if (basehaz == "bs") {
 
    bknots <- c(min_t, max_t)
    iknots <- get_iknots(tt, df = df, iknots = knots, degree = degree)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, degree = degree, type = "bs")      
    nvars  <- ncol(basis)  # number of aux parameters, basis terms
    
  } else if (basehaz == "ms") {
    
    bknots <- c(min_t, max_t)
    iknots <- get_iknots(tt, df = df, iknots = knots, degree = degree, intercept = TRUE)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, degree = degree, type = "ms")      
    nvars  <- ncol(basis)  # number of aux parameters, basis terms
    
  } else if (basehaz == "piecewise") {
    
    degree <- NULL               # degree for splines
    bknots <- c(min_t, max_t)
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- NULL               # spline basis
    nvars  <- length(iknots) + 1 # number of aux parameters, dummy indicators
    
  }
  
  nlist(type_name = basehaz, 
        type = basehaz_for_stan(basehaz),
        nvars, 
        iknots, 
        bknots,
        degree,
        basis,
        df = nvars,
        user_df = nvars,
        knots = if (basehaz == "bs") iknots else c(bknots[1], iknots, bknots[2]),
        bs_basis = basis)
}

# Return a vector with valid names for elements in the list passed to the
# 'basehaz_ops' argument of a 'stan_jm' or 'stan_surv' call
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return A character vector, or NA if unmatched.
get_ok_basehaz_ops <- function(basehaz_name) {
  switch(basehaz_name,
         "bs"        = c("df", "knots", "degree"),
         "ms"        = c("df", "knots", "degree"),
         "piecewise" = c("df", "knots"),
         NA)
}

# Return the integer representation for the baseline hazard, used by Stan
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return An integer, or NA if unmatched.
basehaz_for_stan <- function(basehaz_name) {
  switch(basehaz_name, 
         "weibull"     = 1L,
         "bs"          = 2L,
         "piecewise"   = 3L,
         "ms"          = 4L,
         "exp"         = 5L,
         "gompertz"    = 6L,
         "exp-aft"     = 7L,
         "weibull-aft" = 8L,
         NA)
}

# Return a vector with internal knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - intercept'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is 
#   specified then 'df' is ignored.
# @param degree Non-negative integer. The degree for the spline basis.
# @param iknots Optional vector of internal knots.
# @return A numeric vector of internal knot locations, or NULL if there are
#   no internal knots.
get_iknots <- function(x, df = 5L, degree = 3L, iknots = NULL, intercept = FALSE) {
  
  # obtain number of internal knots
  if (is.null(iknots)) {
    nk <- df - degree - intercept
  } else {
    nk <- length(iknots)
  }
  
  # validate number of internal knots
  if (nk < 0) {
    stop2("Number of internal knots cannot be negative.")
  }
  
  # if no internal knots then return empty vector
  if (nk == 0) {
    return(numeric(0))
  }
  
  # obtain default knot locations if necessary
  if (is.null(iknots)) {
    iknots <- qtile(x, nq = nk + 1) # evenly spaced percentiles
  }
  
  # return internal knot locations, ensuring they are positive
  validate_positive_scalar(iknots)
  
  return(iknots)
}

# Identify whether the type of baseline hazard requires an intercept in
# the linear predictor (NB splines incorporate the intercept into the basis).
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A Logical.
has_intercept <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  (nm %in% c("exp", 
             "exp-aft", 
             "weibull", 
             "weibull-aft", 
             "gompertz", 
             "ms", 
             "bs"))
}

# Return the name of the tve spline coefs or smoothing parameters.
#
# @param x The predictor matrix for the time-varying effects, with column names.
# @param type The type of information about the smoothing parameters to return.
# @return A character or numeric vector, depending on 'type'.
get_smooth_name <- function(x, type = "smooth_coefs") {
  
  if (is.null(x) || !ncol(x))
    return(NULL)  

  nms <- colnames(x)
  nms <- gsub(":splines2::bSpline\\(times__.*\\)[0-9]*$", ":tve-bs-coef", nms)

  nms_trim <- gsub(":tve-[a-z][a-z]-coef[0-9]*$", "", nms)
  tally    <- table(nms_trim)
  indices  <- uapply(tally, seq_len)

  switch(type,
         "smooth_coefs" = paste0(nms, indices),
         "smooth_sd"    = paste0("smooth_sd[", unique(nms_trim), "]"),
         "smooth_map"   = aa(rep(seq_along(tally), tally)),
         "smooth_vars"  = unique(nms_trim),
         stop2("Bug found: invalid input to 'type' argument."))
}

# Return the valid prior distributions for 'prior_aux'.
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A named list.
get_ok_priors_for_aux <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  switch(nm,
         "exp"         = nlist(),
         "exp-aft"     = nlist(),
         "weibull"     = nlist("normal", student_t = "t", "cauchy", "exponential"),
         "weibull-aft" = nlist("normal", student_t = "t", "cauchy", "exponential"),
         "gompertz"    = nlist("normal", student_t = "t", "cauchy", "exponential"),
         "ms"          = nlist("dirichlet"),
         "bs"          = nlist("normal", student_t = "t", "cauchy"),
         "piecewise"   = nlist("normal", student_t = "t", "cauchy"),
         stop2("Bug found: unknown type of baseline hazard."))
}

# Return the default prior distribution for 'prior_aux'.
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A list corresponding to the default prior.
get_default_prior_for_aux <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  switch(nm,
         "exp"         = list(), # equivalent to NULL
         "exp-aft"     = list(), # equivalent to NULL
         "weibull"     = normal(),
         "weibull-aft" = normal(),
         "gompertz"    = normal(),
         "ms"          = dirichlet(),
         "bs"          = normal(),
         "piecewise"   = normal(),
         stop2("Bug found: unknown type of baseline hazard."))
}

# Return the names for the group-specific parameters
#
# @param group List returned by rstanarm:::pad_reTerms.
# @return A character vector.
get_b_names <- function(group) {
  if (is.null(group))
    return(NULL) # no random effects structure
  c(paste0("b[", make_b_nms(group), "]"))
}

# Return the names for the var-cov parameters
#
# @param group List returned by rstanarm:::pad_reTerms.
# @return A character vector.
get_varcov_names <- function(group) {
  if (is.null(group))
    return(NULL) # no random effects structure
  paste0("Sigma[", get_Sigma_nms(group$cnms), "]")
}

# Return the default scale parameter for 'prior_aux'.
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A scalar.
get_default_aux_scale <- function(basehaz) {
  switch(get_basehaz_name(basehaz),
         "weibull"     = 2,
         "weibull-aft" = 2,
         "gompertz"    = 0.5,
         20)
}

# Check if the type of baseline hazard has a closed form
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A logical.
check_for_closed_form <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  nm %in% c("exp",
            "exp-aft",
            "weibull",
            "weibull-aft",
            "gompertz",
            "ms")
}

# Replace the parameter names slot of an object of class 'stanfit'.
#
# @param stanfit An object of class 'stanfit'.
# @param new_nms A character vector of new parameter names.
# @return A 'stanfit' object.
replace_stanfit_nms <- function(stanfit, new_nms) {
  stanfit@sim$fnames_oi <- new_nms
  stanfit
}

# Return the spline basis for the given type of baseline hazard.
# 
# @param times A numeric vector of times at which to evaluate the basis.
# @param basehaz A list with info about the baseline hazard, returned by a 
#   call to 'handle_basehaz'.
# @param integrate A logical, specifying whether to calculate the integral of
#   the specified basis.
# @return A matrix.
make_basis <- function(times, basehaz, integrate = FALSE) {
  N <- length(times)
  K <- basehaz$nvars
  if (!N) { # times is NULL or empty vector
    return(matrix(0, 0, K))
  } 
  switch(basehaz$type_name,
         "exp"         = matrix(0, N, K), # dud matrix for Stan
         "exp-aft"     = matrix(0, N, K), # dud matrix for Stan
         "weibull"     = matrix(0, N, K), # dud matrix for Stan
         "weibull-aft" = matrix(0, N, K), # dud matrix for Stan
         "gompertz"    = matrix(0, N, K), # dud matrix for Stan
         "ms"          = basis_matrix(times, basis = basehaz$basis, integrate = integrate),
         "bs"          = basis_matrix(times, basis = basehaz$basis),
         "piecewise"   = dummy_matrix(times, knots = basehaz$knots),
         stop2("Bug found: unknown type of baseline hazard."))
}

# Evaluate a spline basis matrix at the specified times
#
# @param time A numeric vector.
# @param basis Info on the spline basis.
# @param integrate A logical, should the integral of the basis be returned?
# @return A two-dimensional array.
basis_matrix <- function(times, basis, integrate = FALSE) {
  out <- predict(basis, times)
  if (integrate) {
    stopifnot(inherits(basis, "mSpline"))
    class(basis) <- c("matrix", "iSpline")
    out <- predict(basis, times)
  }
  aa(out)
}

# Parse the model formula and data
#
# @param formula The user input to the formula argument.
# @param data The user input to the data argument (i.e. a data frame).
# @param A list with the model data (following removal of NA rows etc) and
#   a number of elements corresponding to different parts of the formula.
parse_formula_and_data <- function(formula, data) {
  
  formula <- validate_formula(formula, needs_response = TRUE)
  
  # all variables of entire formula
  allvars <- all.vars(formula)
  allvars_form <- reformulate(allvars)
  
  # LHS of entire formula
  lhs       <- lhs(formula)         # LHS as expression
  lhs_form  <- reformulate_lhs(lhs) # LHS as formula
  
  # RHS of entire formula
  rhs       <- rhs(formula)         # RHS as expression
  rhs_form  <- reformulate_rhs(rhs) # RHS as formula

  # evaluate model data (row subsetting etc)
  data <- make_model_data(allvars_form, data)
  
  # evaluated response variables
  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")
  
  if (type == "right") {
    min_t    <- 0
    max_t    <- max(surv[, "time"])
    status   <- as.vector(surv[, "status"])
    t_end    <- as.vector(surv[, "time"])
  } else if (type == "counting") {
    min_t    <- min(surv[, "start"])
    max_t    <- max(surv[, "stop"])
    status   <- as.vector(surv[, "status"])
    t_end    <- as.vector(surv[, "stop"])
  } else if (type == "interval") {
    min_t    <- 0
    max_t    <- max(surv[, c("time1", "time2")])
    status   <- as.vector(surv[, "status"])
    t_end    <- as.vector(surv[, "time1"])
  } else if (type == "interval2") {
    min_t    <- 0
    max_t    <- max(surv[, c("time1", "time2")])
    status   <- as.vector(surv[, "status"])
    t_end    <- as.vector(surv[, "time1"])
  }

  if (any(is.na(status)))
    stop2("Invalid status indicator in Surv object.")
  
  if (any(status < 0 || status > 3))
    stop2("Invalid status indicator in Surv object.")
  
  # deal with tve(x, ...)
  tve_stuff <- handle_tve(formula, 
                          min_t  = min_t, 
                          max_t  = max_t, 
                          times  = t_end, 
                          status = status)
  tf_form    <- tve_stuff$tf_form
  td_form    <- tve_stuff$td_form    # may be NULL
  tt_vars    <- tve_stuff$tt_vars    # may be NULL
  tt_frame   <- tve_stuff$tt_frame   # may be NULL
  tt_types   <- tve_stuff$tt_types   # may be NULL
  tt_degrees <- tve_stuff$tt_degrees # may be NULL
  tt_calls   <- tve_stuff$tt_calls   # may be NULL
  tt_forms   <- tve_stuff$tt_forms   # may be NULL

  # just fixed-effect part of formula
  fe_form   <- lme4::nobars(tf_form)

  # just random-effect part of formula
  bars      <- lme4::findbars(tf_form)
  re_parts  <- lapply(bars, split_at_bars)
  re_forms  <- fetch(re_parts, "re_form")  
  if (length(bars) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")

  nlist(formula,
        data,
        allvars,
        allvars_form,
        lhs,
        lhs_form,
        rhs,
        rhs_form,
        tf_form,
        td_form,
        tt_vars,
        tt_frame,
        tt_types,
        tt_degrees,
        tt_calls,
        tt_forms,
        fe_form,
        bars,
        re_parts,
        re_forms,
        surv_type = attr(surv, "type"))
}

# Handle the 'tve(x, ...)' terms in the model formula
#
# @param Terms terms object for the fixed effect part of the model formula.
# @return A named list with the following elements:
# 
handle_tve <- function(formula, min_t, max_t, times, status) {

  # extract terms objects for fixed effect part of model formula
  Terms <- delete.response(terms(lme4::nobars(formula), specials = "tve"))

  # check which fixed effect terms have a tve() wrapper
  sel <- attr(Terms, "specials")$tve
  
  # if no tve() terms then just return the fixed effect formula as is
  if (!length(sel)) {
    return(list(tf_form  = formula,
                td_form  = NULL,
                tt_vars  = NULL,
                tt_frame = NULL,
                tt_calls = NULL,
                tt_forms = NULL))
  }
   
  # otherwise extract rhs of formula
  all_vars <- rownames(attr(Terms, "factors")) # all variables in fe formula
  tve_vars <- all_vars[sel]                    # variables with a tve() wrapper
  
  # replace 'tve(x, ...)' in formula with 'x'
  old_vars <- all_vars
  new_vars <- sapply(old_vars, function(x) {
    if (x %in% tve_vars) {
      # strip tve() from variable
      tve <- function(y, ...) { safe_deparse(substitute(y)) } # define locally
      return(eval(parse(text = x)))
    } else {
      # just return variable
      return(x)
    }
  }, USE.NAMES = FALSE)
  tf_terms <- attr(Terms, "term.labels")
  td_terms <- c()
  k <- 0 # initialise td_terms indexing (for creating a new formula)
  for (i in sel) {
    sel_terms <- which(attr(Terms, "factors")[i, ] > 0)
    for (j in sel_terms) {
      k <- k + 1
      tf_terms[j] <- td_terms[k] <- gsub(old_vars[i], 
                                         new_vars[i], 
                                         tf_terms[j], 
                                         fixed = TRUE)
    }
  }
  
  # extract 'tve(x, ...)' from formula and return '~ x' and '~ bs(times, ...)'
  idx <- 1
  tt_vars    <- list()
  tt_types   <- list()
  tt_degrees <- list()
  tt_calls   <- list()
  
  for (i in seq_along(sel)) {
    
    # define tve() function locally; uses objects from the parent environment
    #
    # @param x The variable the time-varying effect is going to be applied to.
    # @param type Character string, the type of time-varying effect to use. Can
    #   currently only be "bs".
    # @param df,knots,degree Additional arguments passed to splines2::bSpline.
    # @return The call used to construct a time-varying basis.
    tve <- function(x, type = "bs", df = NULL, knots = NULL, degree = 3L) {
      
      type <- match.arg(type)
      
      if (!is.null(df) && !is.null(knots))
        stop("Cannot specify both 'df' and 'knots' in the 'tve' function.")
      
      if (degree < 0)
        stop("In 'tve' function, 'degree' must be non-negative.")
      
      if (is.null(df) && is.null(knots))
        df <- 3L
            
      # note that times and status are taken from the parent environment
      tt <- times[status == 1] # uncensored event times
      if (is.null(knots) && !length(tt)) {
        warning2("No observed events found in the data. Censoring times will ",
                 "be used to evaluate default knot locations for tve().")
        tt <- times
      }
      
      # note that min_t and max_t are taken from the parent environment
      if (!is.null(knots)) {
        if (any(knots < min_t))
          stop2("In tve(), 'knots' cannot be placed before the earliest entry time.")
        if (any(knots > max_t))
          stop2("In tve(), 'knots' cannot be placed beyond the latest event time.")
      }

      if (type == "bs") {
        
        iknots <- get_iknots(tt, df = df, iknots = knots)
 
        bknots <- c(min_t, max_t)
        
        new_args <- list(knots          = iknots,
                         Boundary.knots = bknots,
                         degree         = degree)

        return(list(
          type   = type,
          degree = degree,
          call   = sub("^list\\(", "splines2::bSpline\\(times__, ", 
                       deparse(new_args, 500L, control = c("all", "hexNumeric")))))
                       # NB use of hexNumeric to ensure numeric accuracy is maintained
        
      }

    }
    
    tt_parsed <- eval(parse(text = all_vars[sel[i]]))
    tt_terms  <- which(attr(Terms, "factors")[i, ] > 0)
    for (j in tt_terms) {
      tt_vars   [[idx]] <- tf_terms[j]
      tt_types  [[idx]] <- tt_parsed$type
      tt_degrees[[idx]] <- tt_parsed$degree
      tt_calls  [[idx]] <- tt_parsed$call
      idx <- idx + 1
    }
  }
  
  # add on the terms labels from the random effects part of the formula
  bars <- lme4::findbars(formula)
  if (length(bars)) {
    re_terms <- sapply(bars, bracket_wrap)
    tf_terms <- c(tf_terms, re_terms)
  }
  
  # formula with all variables but no 'tve(x, ...)' wrappers
  tf_form <- reformulate(tf_terms, response = lhs(formula))
  
  # formula with only tve variables but no 'tve(x, ...)' wrappers
  td_form <- reformulate(td_terms, response = lhs(formula))
  
  # unique set of '~ bs(times__, ...)' calls based on all 'tve(x, ...)' terms
  tt_frame <- reformulate(unique(unlist(tt_calls)), intercept = FALSE)
  
  # formula with '~ x' and '~ bs(times__, ...)' from each 'tve(x, ...)' call
  tt_vars  <- lapply(tt_vars,  reformulate)
  tt_forms <- lapply(tt_calls, reformulate)
  
  # return object
  nlist(tf_form,
        td_form,
        tt_vars,
        tt_frame,
        tt_types,
        tt_degrees,
        tt_calls,
        tt_forms)
}

# Ensure only valid arguments are passed to the tve() call
validate_tve_args <- function(dots, ok_args) {
  
  if (!isTRUE(all(names(dots) %in% ok_args)))
    stop2("Invalid argument to 'tve' function. ",
          "Valid arguments are: ", comma(ok_args))  
  
  return(dots)
}

# Deparse an expression and wrap it in brackets
#
# @param x An expression.
# @return A character string.
bracket_wrap <- function(x) {
  paste0("(", deparse(x, 500), ")")
}

# Check input to the formula argument
#
# @param formula The user input to the formula argument.
# @param needs_response A logical; if TRUE then formula must contain a LHS.
# @return A formula.
validate_formula <- function(formula, needs_response = TRUE) {
  
  if (!inherits(formula, "formula")) {
    stop2("'formula' must be a formula.")
  }
  
  if (needs_response) {
    len <- length(formula)
    if (len < 3) {
      stop2("'formula' must contain a response.")
    }
  }
  as.formula(formula)
}

# Check object is a Surv object with a valid type
#
# @param x A Surv object. That is, the LHS of a formula as evaluated in a 
#   data frame environment.
# @param ok_types A character vector giving the valid types of Surv object.
# @return A Surv object.
validate_surv <- function(x, ok_types = c("right", "counting",
                                          "interval", "interval2")) {
  if (!inherits(x, "Surv"))
    stop2("LHS of 'formula' must be a 'Surv' object.")
  if (!attr(x, "type") %in% ok_types)
    stop2("Surv object type must be one of: ", comma(ok_types))
  x
}

# Extract LHS of a formula
#
# @param x A formula object.
# @return An expression.
lhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[2L]]
  } else {
    out <- NULL
  }
  out
}

# Extract RHS of a formula
#
# @param x A formula object.
# @return An expression.
rhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[3L]]
  } else {
    out <- x[[2L]]
  }
  out
}

# Reformulate as LHS of a formula
#
# @param x A character string or expression.
# @return A formula.
reformulate_lhs <- function(x) {
  x <- formula(substitute(LHS ~ 1, list(LHS = x)))
  x
}

# Reformulate as RHS of a formula
#
# @param x A character string or expression.
# @return A formula.
reformulate_rhs <- function(x) {
  x <- formula(substitute(~ RHS, list(RHS = x)))
  x
}

# Return the response vector (time)
#
# @param model_frame The model frame.
# @param type The type of time variable to return:
#   "beg": the entry time for the row in the survival data,
#   "end": the exit  time for the row in the survival data,
#   "gap": the difference between entry and exit times,
#   "upp": if the row involved interval censoring, then the exit time
#          would have been the lower limit of the interval, and "upp" 
#          is the upper limit of the interval.
# @return A numeric vector.
make_t <- function(model_frame, type = c("beg", "end", "gap", "upp")) {
  
  type <- match.arg(type)
  resp <- if (survival::is.Surv(model_frame)) 
    model_frame else model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")
  
  t_beg <- switch(surv,
                  "right"     = rep(0, nrow(model_frame)),
                  "interval"  = rep(0, nrow(model_frame)),
                  "interval2" = rep(0, nrow(model_frame)),
                  "counting"  = as.vector(resp[, "start"]),
                  stop(err))

  t_end <- switch(surv,
                  "right"     = as.vector(resp[, "time"]),
                  "interval"  = as.vector(resp[, "time1"]),
                  "interval2" = as.vector(resp[, "time1"]),
                  "counting"  = as.vector(resp[, "stop"]),
                  stop(err))

  t_upp <- switch(surv,
                  "right"     = rep(NaN, nrow(model_frame)),
                  "counting"  = rep(NaN, nrow(model_frame)),
                  "interval"  = as.vector(resp[, "time2"]),
                  "interval2" = as.vector(resp[, "time2"]),
                  stop(err))

  switch(type,
         "beg" = t_beg,
         "end" = t_end,
         "gap" = t_end - t_beg,
         "upp" = t_upp,
         stop("Bug found: cannot handle specified 'type'."))
}

# Return the response vector (status indicator)
#
# @param model_frame The model frame.
# @return A numeric vector.
make_d <- function(model_frame) {
  
  resp <- if (survival::is.Surv(model_frame)) 
    model_frame else model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")
  
  switch(surv,
         "right"     = as.vector(resp[, "status"]),
         "interval"  = as.vector(resp[, "status"]),
         "interval2" = as.vector(resp[, "status"]),
         "counting"  = as.vector(resp[, "status"]),
         stop(err))
}

# Return a data frame with NAs excluded
#
# @param formula The parsed model formula.
# @param data The (user-specified) data frame.
# @return A data frame, with only complete cases for the variables that
#   appear in the model formula.
make_model_data <- function(formula, data) {
  mf <- model.frame(formula, data, na.action = na.pass)
  include <- apply(mf, 1L, function(row) !any(is.na(row)))
  data[include, , drop = FALSE]
}

# Return the model frame
#
# @param formula The parsed model formula.
# @param data The model data frame.
# @param xlevs Passed to xlev argument of model.frame.
# @param drop.unused.levels Passed to drop.unused.levels argument of model.frame.
# @param check_constant If TRUE then an error is thrown is the returned
#   model frame contains any constant variables.
# @return A list with the following elements:
#   mf: the model frame based on the formula.
#   mt: the model terms associated with the returned model frame.
make_model_frame <- function(formula, 
                             data, 
                             xlevs              = NULL,
                             drop.unused.levels = FALSE,
                             check_constant     = FALSE,
                             na.action          = na.fail) {

  # construct model frame
  Terms <- terms(lme4::subbars(formula))
  mf <- stats::model.frame(Terms, 
                           data,
                           xlev = xlevs, 
                           drop.unused.levels = drop.unused.levels,
                           na.action = na.action)
  
  # get predvars for fixed part of formula
  TermsF <- terms(lme4::nobars(formula)) 
  mfF <- stats::model.frame(TermsF, 
                            data, 
                            xlev = xlevs, 
                            drop.unused.levels = drop.unused.levels,
                            na.action = na.action)
  attr(attr(mf, "terms"), "predvars.fixed") <- attr(attr(mfF, "terms"), "predvars")
  
  # get predvars for random part of formula
  has_bars <- length(lme4::findbars(formula)) > 0
  if (has_bars) {
    TermsR <- terms(lme4::subbars(justRE(formula, response = TRUE)))
    mfR <- stats::model.frame(TermsR,
                              data, 
                              xlev = xlevs, 
                              drop.unused.levels = drop.unused.levels,
                              na.action = na.action)
    attr(attr(mf, "terms"), "predvars.random") <- attr(attr(mfR, "terms"), "predvars")
  } else {
    attr(attr(mf, "terms"), "predvars.random") <- NULL
  }
  
  # check no constant vars
  if (check_constant)
    mf <- check_constant_vars(mf)
  
  # add additional predvars attributes
  
  # check for terms
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) 
    stop2("No intercept or predictors specified.")
  
  nlist(mf, mt)
}

# Return the predictor matrix
#
# @param formula The parsed model formula.
# @param model_frame The model frame.
# @param xlevs Passed to xlev argument of model.matrix.
# @param check_constant If TRUE then an error is thrown is the returned
#   predictor matrix contains any constant columns.
# @return A named list with the following elements:
#   x: the fe model matrix, not centered and without intercept.
#   x_bar: the column means of the model matrix.
#   x_centered: the fe model matrix, centered.
#   N: number of rows (observations) in the model matrix.
#   K: number of cols (predictors) in the model matrix.
make_x <- function(formula, 
                   model_frame, 
                   xlevs = NULL,
                   check_constant = TRUE) {

  # uncentred predictor matrix, without intercept
  x <- model.matrix(formula, model_frame, xlev = xlevs)
  x <- drop_intercept(x)
  
  # column means of predictor matrix
  x_bar <- aa(colMeans(x))

  # centered predictor matrix
  x_centered <- sweep(x, 2, x_bar, FUN = "-")
  
  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (apply(x, 2L, n_distinct) < 2)
  if (check_constant && any(sel)) {
    cols <- paste(colnames(x)[sel], collapse = ", ")
    stop2("Cannot deal with empty interaction levels found in columns: ", cols)
  }
  
  nlist(x, x_centered, x_bar, N = NROW(x), K = NCOL(x))
}

# Return the tve predictor matrix
#
# @param formula The parsed model formula.
# @param model_frame The model frame.
# @param xlevs Passed to xlev argument of model.matrix.
# @return A named list with the following elements:
#   s: model matrix for time-varying terms, not centered and without intercept.
#   tt_ncol: stored attribute, a numeric vector with the number of columns in
#     the model matrix that correspond to each tve() term in the original 
#     model formula.
#   tt_map: stored attribute, a numeric vector with indexing for the columns 
#     of the model matrix stating which tve() term in the original model 
#     formula they correspond to.
make_s <- function(formula,
                   model_frame, 
                   xlevs = NULL) {
  
  # create the design matrix for each tve() term
  s_parts <- xapply(
    formula$tt_vars,  # names of variables with a tve() wrapper
    formula$tt_forms, # time-transformation functions to interact them with
    FUN = function(vn, tt) {
      m1 <- make_x(vn, model_frame, xlevs = xlevs, check_constant = FALSE)$x
      m2 <- make_x(tt, model_frame, xlevs = xlevs, check_constant = FALSE)$x
      m3 <- matrix(apply(m1, 2L, `*`, m2), nrow = nrow(m2))
      colnames(m3) <- uapply(colnames(m1), paste, colnames(m2), sep = ":")
      return(m3)
    })
  
  # bind columns to form one design matrix for tve() terms
  s <- do.call("cbind", s_parts)
  
  # store indexing of the columns in the design matrix
  tt_ncol <- sapply(s_parts, ncol)
  tt_map  <- rep(seq_along(tt_ncol), tt_ncol)
  
  # return design matrix with indexing info as an attribute
  structure(s, tt_ncol = tt_ncol, tt_map  = tt_map)
}

# Check if the only element of a character vector is 'Intercept'
#
# @param x A character vector.
# @return A logical.
intercept_only <- function(x) {
  length(x) == 1 && x == "(Intercept)"
}
