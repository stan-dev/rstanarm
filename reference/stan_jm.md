# Bayesian joint longitudinal and time-to-event models via Stan

Fits a shared parameter joint model for longitudinal and time-to-event
(e.g. survival) data under a Bayesian framework using Stan.

## Usage

``` r
stan_jm(
  formulaLong,
  dataLong,
  formulaEvent,
  dataEvent,
  time_var,
  id_var,
  family = gaussian,
  assoc = "etavalue",
  lag_assoc = 0,
  grp_assoc,
  scale_assoc = NULL,
  epsilon = 1e-05,
  basehaz = c("bs", "weibull", "piecewise"),
  basehaz_ops,
  qnodes = 15,
  init = "prefit",
  weights,
  priorLong = normal(autoscale = TRUE),
  priorLong_intercept = normal(autoscale = TRUE),
  priorLong_aux = cauchy(0, 5, autoscale = TRUE),
  priorEvent = normal(autoscale = TRUE),
  priorEvent_intercept = normal(autoscale = TRUE),
  priorEvent_aux = cauchy(autoscale = TRUE),
  priorEvent_assoc = normal(autoscale = TRUE),
  prior_covariance = lkj(autoscale = TRUE),
  prior_PD = FALSE,
  algorithm = c("sampling", "meanfield", "fullrank"),
  adapt_delta = NULL,
  max_treedepth = 10L,
  QR = FALSE,
  sparse = FALSE,
  ...
)
```

## Arguments

- formulaLong:

  A two-sided linear formula object describing both the fixed-effects
  and random-effects parts of the longitudinal submodel, similar in vein
  to formula specification in the **lme4** package (see
  [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) or the **lme4**
  vignette for details). Note however that the double bar (`||`)
  notation is not allowed when specifying the random-effects parts of
  the formula, and neither are nested grouping factors (e.g.
  `(1 | g1/g2))` or `(1 | g1:g2)`, where `g1`, `g2` are grouping
  factors. Offset terms can also be included in the model formula. For a
  multivariate joint model (i.e. more than one longitudinal marker) this
  should be a list of such formula objects, with each element of the
  list providing the formula for one of the longitudinal submodels.

- dataLong:

  A data frame containing the variables specified in `formulaLong`. If
  fitting a multivariate joint model, then this can be either a single
  data frame which contains the data for all longitudinal submodels, or
  it can be a list of data frames where each element of the list
  provides the data for one of the longitudinal submodels.

- formulaEvent:

  A two-sided formula object describing the event submodel. The left
  hand side of the formula should be a
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object. See
  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html).

- dataEvent:

  A data frame containing the variables specified in `formulaEvent`.

- time_var:

  A character string specifying the name of the variable in `dataLong`
  which represents time.

- id_var:

  A character string specifying the name of the variable in `dataLong`
  which distinguishes between individuals. This can be left unspecified
  if there is only one grouping factor (which is assumed to be the
  individual). If there is more than one grouping factor (i.e.
  clustering beyond the level of the individual) then the `id_var`
  argument must be specified.

- family:

  The family (and possibly also the link function) for the longitudinal
  submodel(s). See [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html)
  for details. If fitting a multivariate joint model, then this can
  optionally be a list of families, in which case each element of the
  list specifies the family for one of the longitudinal submodels.

- assoc:

  A character string or character vector specifying the joint model
  association structure. Possible association structures that can be
  used include: "etavalue" (the default); "etaslope"; "etaauc";
  "muvalue"; "muslope"; "muauc"; "shared_b"; "shared_coef"; or "null".
  These are described in the **Details** section below. For a
  multivariate joint model, different association structures can
  optionally be used for each longitudinal submodel by specifying a list
  of character vectors, with each element of the list specifying the
  desired association structure for one of the longitudinal submodels.
  Specifying `assoc = NULL` will fit a joint model with no association
  structure (equivalent to fitting separate longitudinal and
  time-to-event models). It is also possible to include interaction
  terms between the association term ("etavalue", "etaslope", "muvalue",
  "muslope") and observed data/covariates. It is also possible, when
  fitting a multivariate joint model, to include interaction terms
  between the association terms ("etavalue" or "muvalue") corresponding
  to the different longitudinal outcomes. See the **Details** section as
  well as the **Examples** below.

- lag_assoc:

  A non-negative scalar specifying the time lag that should be used for
  the association structure. That is, the hazard of the event at time
  *t* will be assumed to be associated with the value/slope/auc of the
  longitudinal marker at time *t-u*, where *u* is the time lag. If
  fitting a multivariate joint model, then a different time lag can be
  used for each longitudinal marker by providing a numeric vector of
  lags, otherwise if a scalar is provided then the specified time lag
  will be used for all longitudinal markers. Note however that only one
  time lag can be specified for linking each longitudinal marker to the
  event, and that that time lag will be used for all association
  structure types (e.g. `"etavalue"`, `"etaslope"`, `"etaauc"`,
  `"muvalue"`, etc) that are specified for that longitudinal marker in
  the `assoc` argument.

- grp_assoc:

  Character string specifying the method for combining information
  across lower level units clustered within an individual when forming
  the association structure. This is only relevant when a grouping
  factor is specified in `formulaLong` that corresponds to clustering
  within individuals. This can be specified as either `"sum"`, `mean`,
  `"min"` or `"max"`. For example, specifying `grp_assoc = "sum"`
  indicates that the association structure should be based on a
  summation across the lower level units clustered within an individual,
  or specifying `grp_assoc = "mean"` indicates that the association
  structure should be based on the mean (i.e. average) taken across the
  lower level units clustered within an individual. So, for example,
  specifying `assoc = "muvalue"` and `grp_assoc = "sum"` would mean that
  the log hazard at time *t* for individual *i* would be linearly
  related to the sum of the expected values at time *t* for each of the
  lower level units (which may be for example tumor lesions) clustered
  within that individual.

- scale_assoc:

  A non-zero numeric value specifying an optional scaling parameter for
  the association structure. This multiplicatively scales the
  value/slope/auc of the longitudinal marker by `scale_assoc` within the
  event submodel. When fitting a multivariate joint model, a scaling
  parameter must be specified for each longitudinal submodel using a
  vector of numeric values. Note that only one scaling parameter can be
  specified for each longitudinal submodel, and it will be used for all
  association structure types (e.g. `"etavalue"`, `"etaslope"`,
  `"etaauc"`, `"muvalue"`, etc) that are specified for that longitudinal
  marker in the `assoc` argument.

- epsilon:

  The half-width of the central difference used to numerically calculate
  the derivate when the `"etaslope"` association structure is used.

- basehaz:

  A character string indicating which baseline hazard to use for the
  event submodel. Options are a B-splines approximation estimated for
  the log baseline hazard (`"bs"`, the default), a Weibull baseline
  hazard (`"weibull"`), or a piecewise constant baseline hazard
  (`"piecewise"`). (Note however that there is currently limited
  post-estimation functionality available for models estimated using a
  piecewise constant baseline hazard).

- basehaz_ops:

  A named list specifying options related to the baseline hazard.
  Currently this can include:  

  `df`

  :   A positive integer specifying the degrees of freedom for the
      B-splines if `basehaz = "bs"`, or the number of intervals used for
      the piecewise constant baseline hazard if `basehaz = "piecewise"`.
      The default is 6.

  `knots`

  :   An optional numeric vector specifying the internal knot locations
      for the B-splines if `basehaz = "bs"`, or the internal cut-points
      for defining intervals of the piecewise constant baseline hazard
      if `basehaz = "piecewise"`. Knots cannot be specified if `df` is
      specified. If not specified, then the default is to use `df - 4`
      knots if `basehaz = "bs"`, or `df - 1` knots if
      `basehaz = "piecewise"`, which are placed at equally spaced
      percentiles of the distribution of observed event times.

- qnodes:

  The number of nodes to use for the Gauss-Kronrod quadrature that is
  used to evaluate the cumulative hazard in the likelihood function.
  Options are 15 (the default), 11 or 7.

- init:

  The method for generating the initial values for the MCMC. The default
  is `"prefit"`, which uses those obtained from fitting separate
  longitudinal and time-to-event models prior to fitting the joint
  model. The separate longitudinal model is a (possibly multivariate)
  generalised linear mixed model estimated using variational bayes. This
  is achieved via the
  [`stan_mvmer`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md)
  function with `algorithm = "meanfield"`. The separate Cox model is
  estimated using
  [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html). This is
  achieved using the and time-to-event models prior to fitting the joint
  model. The separate models are estimated using the
  [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) and
  [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) functions. This
  should provide reasonable initial values which should aid the MCMC
  sampler. Parameters that cannot be obtained from fitting separate
  longitudinal and time-to-event models are initialised using the
  "random" method for
  [`stan`](https://mc-stan.org/rstan/reference/stan.html). However it is
  recommended that any final analysis should ideally be performed with
  several MCMC chains each initiated from a different set of initial
  values; this can be obtained by setting `init = "random"`. In
  addition, other possibilities for specifying `init` are the same as
  those described for
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- weights:

  Experimental and should be used with caution. The user can optionally
  supply a 2-column data frame containing a set of 'prior weights' to be
  used in the estimation process. The data frame should contain two
  columns: the first containing the IDs for each individual, and the
  second containing the corresponding weights. The data frame should
  only have one row for each individual; that is, weights should be
  constant within individuals.

- priorLong, priorEvent, priorEvent_assoc:

  The prior distributions for the regression coefficients in the
  longitudinal submodel(s), event submodel, and the association
  parameter(s). Can be a call to one of the various functions provided
  by rstanarm for specifying priors. The subset of these functions that
  can be used for the prior on the coefficients can be grouped into
  several "families":

  |                                 |                                 |
  |---------------------------------|---------------------------------|
  | **Family**                      | **Functions**                   |
  | *Student t family*              | `normal`, `student_t`, `cauchy` |
  | *Hierarchical shrinkage family* | `hs`, `hs_plus`                 |
  | *Laplace family*                | `laplace`, `lasso`              |

  See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) for details on
  the families and how to specify the arguments for all of the functions
  in the table above. To omit a prior —i.e., to use a flat (improper)
  uniform prior— `prior` can be set to `NULL`, although this is rarely a
  good idea.

  **Note:** Unless `QR=TRUE`, if `prior` is from the Student t family or
  Laplace family, and if the `autoscale` argument to the function used
  to specify the prior (e.g.
  [`normal`](https://mc-stan.org/rstanarm/reference/priors.md)) is left
  at its default and recommended value of `TRUE`, then the default or
  user-specified prior scale(s) may be adjusted internally based on the
  scales of the predictors. See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) for details on
  the rescaling and the
  [`prior_summary`](https://mc-stan.org/rstanarm/reference/prior_summary.stanreg.md)
  function for a summary of the priors used for a particular model.

- priorLong_intercept, priorEvent_intercept:

  The prior distributions for the intercepts in the longitudinal
  submodel(s) and event submodel. Can be a call to `normal`, `student_t`
  or `cauchy`. See the [priors help
  page](https://mc-stan.org/rstanarm/reference/priors.md) for details on
  these functions. To omit a prior on the intercept —i.e., to use a flat
  (improper) uniform prior— `prior_intercept` can be set to `NULL`.

  **Note:** The prior distribution for the intercept is set so it
  applies to the value when all predictors are centered. Moreover, note
  that a prior is only placed on the intercept for the event submodel
  when a Weibull baseline hazard has been specified. For the B-splines
  and piecewise constant baseline hazards there is not intercept
  parameter that is given a prior distribution; an intercept parameter
  will be shown in the output for the fitted model, but this just
  corresponds to the necessary post-estimation adjustment in the linear
  predictor due to the centering of the predictiors in the event
  submodel.

- priorLong_aux:

  The prior distribution for the "auxiliary" parameters in the
  longitudinal submodels (if applicable). The "auxiliary" parameter
  refers to a different parameter depending on the `family`. For
  Gaussian models `priorLong_aux` controls `"sigma"`, the error standard
  deviation. For negative binomial models `priorLong_aux` controls
  `"reciprocal_dispersion"`, which is similar to the `"size"` parameter
  of [`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html): smaller
  values of `"reciprocal_dispersion"` correspond to greater dispersion.
  For gamma models `priorLong_aux` sets the prior on to the `"shape"`
  parameter (see e.g.,
  [`rgamma`](https://rdrr.io/r/stats/GammaDist.html)), and for
  inverse-Gaussian models it is the so-called `"lambda"` parameter
  (which is essentially the reciprocal of a scale parameter). Binomial
  and Poisson models do not have auxiliary parameters.

  `priorLong_aux` can be a call to `exponential` to use an exponential
  distribution, or `normal`, `student_t` or `cauchy`, which results in a
  half-normal, half-t, or half-Cauchy prior. See
  [`priors`](https://mc-stan.org/rstanarm/reference/priors.md) for
  details on these functions. To omit a prior —i.e., to use a flat
  (improper) uniform prior— set `priorLong_aux` to `NULL`.

  If fitting a multivariate joint model, you have the option to specify
  a list of prior distributions, however the elements of the list that
  correspond to any longitudinal submodel which does not have an
  auxiliary parameter will be ignored.

- priorEvent_aux:

  The prior distribution for the "auxiliary" parameters in the event
  submodel. The "auxiliary" parameters refers to different parameters
  depending on the baseline hazard. For `basehaz = "weibull"` the
  auxiliary parameter is the Weibull shape parameter. For
  `basehaz = "bs"` the auxiliary parameters are the coefficients for the
  B-spline approximation to the log baseline hazard. For
  `basehaz = "piecewise"` the auxiliary parameters are the piecewise
  estimates of the log baseline hazard.

- prior_covariance:

  Cannot be `NULL`; see
  [`priors`](https://mc-stan.org/rstanarm/reference/priors.md) for more
  information about the prior distributions on covariance matrices. Note
  however that the default prior for covariance matrices in `stan_jm` is
  slightly different to that in
  [`stan_glmer`](https://mc-stan.org/rstanarm/reference/stan_glmer.md)
  (the details of which are described on the
  [`priors`](https://mc-stan.org/rstanarm/reference/priors.md) page).

- prior_PD:

  A logical scalar (defaulting to `FALSE`) indicating whether to draw
  from the prior predictive distribution instead of conditioning on the
  outcome.

- algorithm:

  A string (possibly abbreviated) indicating the estimation approach to
  use. Can be `"sampling"` for MCMC (the default), `"optimizing"` for
  optimization, `"meanfield"` for variational inference with independent
  normal distributions, or `"fullrank"` for variational inference with a
  multivariate normal distribution. See
  [`rstanarm-package`](https://mc-stan.org/rstanarm/reference/rstanarm-package.md)
  for more details on the estimation algorithms. NOTE: not all fitting
  functions support all four algorithms.

- adapt_delta:

  Only relevant if `algorithm="sampling"`. See the
  [adapt_delta](https://mc-stan.org/rstanarm/reference/adapt_delta.md)
  help page for details.

- max_treedepth:

  A positive integer specifying the maximum treedepth for the non-U-turn
  sampler. See the `control` argument in
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- QR:

  A logical scalar defaulting to `FALSE`, but if `TRUE` applies a scaled
  [`qr`](https://rdrr.io/r/base/qr.html) decomposition to the design
  matrix. The transformation does not change the likelihood of the data
  but is recommended for computational reasons when there are multiple
  predictors. See the
  [QR-argument](https://mc-stan.org/rstanarm/reference/QR-argument.md)
  documentation page for details on how rstanarm does the transformation
  and important information about how to interpret the prior
  distributions of the model parameters when using `QR=TRUE`.

- sparse:

  A logical scalar (defaulting to `FALSE`) indicating whether to use a
  sparse representation of the design (X) matrix. If `TRUE`, the the
  design matrix is not centered (since that would destroy the sparsity)
  and likewise it is not possible to specify both `QR = TRUE` and
  `sparse = TRUE`. Depending on how many zeros there are in the design
  matrix, setting `sparse = TRUE` may make the code run faster and can
  consume much less RAM.

- ...:

  Further arguments passed to the function in the rstan package
  (`sampling`, `vb`, or `optimizing`), corresponding to the estimation
  method named by `algorithm`. For example, if `algorithm` is
  `"sampling"` it is possible to specify `iter`, `chains`, `cores`, and
  other MCMC controls.

  Another useful argument that can be passed to rstan via `...` is
  `refresh`, which specifies how often to print updates when sampling
  (i.e., show the progress every `refresh` iterations). `refresh=0`
  turns off the iteration updates.

## Value

A [stanjm](https://mc-stan.org/rstanarm/reference/stanreg-objects.md)
object is returned.

## Details

The `stan_jm` function can be used to fit a joint model (also known as a
shared parameter model) for longitudinal and time-to-event data under a
Bayesian framework. The underlying estimation is carried out using the
Bayesian C++ package Stan (<https://mc-stan.org/>).  
  
The joint model may be univariate (with only one longitudinal submodel)
or multivariate (with more than one longitudinal submodel). For the
longitudinal submodel a (possibly multivariate) generalised linear mixed
model is assumed with any of the
[`family`](https://rdrr.io/r/stats/family.html) choices allowed by
[`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html). If a multivariate
joint model is specified (by providing a list of formulas in the
`formulaLong` argument), then the multivariate longitudinal submodel
consists of a multivariate generalized linear model (GLM) with
group-specific terms that are assumed to be correlated across the
different GLM submodels. That is, within a grouping factor (for example,
patient ID) the group-specific terms are assumed to be correlated across
the different GLM submodels. It is possible to specify a different
outcome type (for example a different family and/or link function) for
each of the GLM submodels, by providing a list of
[`family`](https://rdrr.io/r/stats/family.html) objects in the `family`
argument. Multi-level clustered data are allowed, and that additional
clustering can occur at a level higher than the individual-level (e.g.
patients clustered within clinics), or at a level lower than the
individual-level (e.g. tumor lesions clustered within patients). If the
clustering occurs at a level lower than the individual, then the user
needs to indicate how the lower level clusters should be handled when
forming the association structure between the longitudinal and event
submodels (see the `grp_assoc` argument described above).  
  
For the event submodel a parametric proportional hazards model is
assumed. The baseline hazard can be estimated using either a cubic
B-splines approximation (`basehaz = "bs"`, the default), a Weibull
distribution (`basehaz = "weibull"`), or a piecewise constant baseline
hazard (`basehaz = "piecewise"`). If the B-spline or piecewise constant
baseline hazards are used, then the degrees of freedom or the internal
knot locations can be (optionally) specified. If the degrees of freedom
are specified (through the `df` argument) then the knot locations are
automatically generated based on the distribution of the observed event
times (not including censoring times). Otherwise internal knot locations
can be specified directly through the `knots` argument. If neither `df`
or `knots` is specified, then the default is to set `df` equal to 6. It
is not possible to specify both `df` and `knots`.  
  
Time-varying covariates are allowed in both the longitudinal and event
submodels. These should be specified in the data in the same way as they
normally would when fitting a separate longitudinal model using
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) or a separate
time-to-event model using
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html). These
time-varying covariates should be exogenous in nature, otherwise they
would perhaps be better specified as an additional outcome (i.e. by
including them as an additional longitudinal outcome in the joint
model).  
  
Bayesian estimation of the joint model is performed via MCMC. The
Bayesian model includes independent priors on the regression
coefficients for both the longitudinal and event submodels, including
the association parameter(s) (in much the same way as the regression
parameters in
[`stan_glm`](https://mc-stan.org/rstanarm/reference/stan_glm.md)) and
priors on the terms of a decomposition of the covariance matrices of the
group-specific parameters. See
[`priors`](https://mc-stan.org/rstanarm/reference/priors.md) for more
information about the priors distributions that are available.  
  
Gauss-Kronrod quadrature is used to numerically evaluate the integral
over the cumulative hazard in the likelihood function for the event
submodel. The accuracy of the numerical approximation can be controlled
using the number of quadrature nodes, specified through the `qnodes`
argument. Using a higher number of quadrature nodes will result in a
more accurate approximation.

### Association structures

The association structure for the joint model can be based on any of the
following parameterisations:

- current value of the linear predictor in the longitudinal submodel
  (`"etavalue"`)

- first derivative (slope) of the linear predictor in the longitudinal
  submodel (`"etaslope"`)

- the area under the curve of the linear predictor in the longitudinal
  submodel (`"etaauc"`)

- current expected value of the longitudinal submodel (`"muvalue"`)

- the area under the curve of the expected value from the longitudinal
  submodel (`"muauc"`)

- shared individual-level random effects (`"shared_b"`)

- shared individual-level random effects which also incorporate the
  corresponding fixed effect as well as any corresponding random effects
  for clustering levels higher than the individual) (`"shared_coef"`)

- interactions between association terms and observed data/covariates
  (`"etavalue_data"`, `"etaslope_data"`, `"muvalue_data"`,
  `"muslope_data"`). These are described further below.

- interactions between association terms corresponding to different
  longitudinal outcomes in a multivariate joint model
  (`"etavalue_etavalue(#)"`, `"etavalue_muvalue(#)"`,
  `"muvalue_etavalue(#)"`, `"muvalue_muvalue(#)"`). These are described
  further below.

- no association structure (equivalent to fitting separate longitudinal
  and event models) (`"null"` or `NULL`)

More than one association structure can be specified, however, not all
possible combinations are allowed. Note that for the lagged association
structures baseline values (time = 0) are used for the instances where
the time lag results in a time prior to baseline. When using the
`"etaauc"` or `"muauc"` association structures, the area under the curve
is evaluated using Gauss-Kronrod quadrature with 15 quadrature nodes. By
default, `"shared_b"` and `"shared_coef"` contribute all random effects
to the association structure; however, a subset of the random effects
can be chosen by specifying their indices between parentheses as a
suffix, for example, `"shared_b(1)"` or `"shared_b(1:3)"` or
`"shared_b(1,2,4)"`, and so on.  
  
In addition, several association terms (`"etavalue"`, `"etaslope"`,
`"muvalue"`, `"muslope"`) can be interacted with observed
data/covariates. To do this, use the association term's main handle plus
a suffix of `"_data"` then followed by the model matrix formula in
parentheses. For example if we had a variable in our dataset for gender
named `sex` then we might want to obtain different estimates for the
association between the current slope of the marker and the risk of the
event for each gender. To do this we would specify
`assoc = c("etaslope", "etaslope_data(~ sex)")`.  
  
It is also possible, when fitting a multivariate joint model, to include
interaction terms between the association terms themselves (this only
applies for interacting `"etavalue"` or `"muvalue"`). For example, if we
had a joint model with two longitudinal markers, we could specify
`assoc = list(c("etavalue", "etavalue_etavalue(2)"), "etavalue")`. The
first element of list says we want to use the value of the linear
predictor for the first marker, as well as it's interaction with the
value of the linear predictor for the second marker. The second element
of the list says we want to also include the expected value of the
second marker (i.e. as a "main effect"). Therefore, the linear predictor
for the event submodel would include the "main effects" for each marker
as well as their interaction.  
  
There are additional examples in the **Examples** section below.

## See also

[`stanreg-objects`](https://mc-stan.org/rstanarm/reference/stanreg-objects.md),
[`stanmvreg-methods`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md),
[`print.stanmvreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md),
[`summary.stanmvreg`](https://mc-stan.org/rstanarm/reference/summary.stanreg.md),
[`posterior_traj`](https://mc-stan.org/rstanarm/reference/posterior_traj.md),
[`posterior_survfit`](https://mc-stan.org/rstanarm/reference/posterior_survfit.md),
[`posterior_predict`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.md),
[`posterior_interval`](https://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.md),
[`pp_check`](https://mc-stan.org/rstanarm/reference/pp_check.stanreg.md),
[`ps_check`](https://mc-stan.org/rstanarm/reference/ps_check.md),
[`stan_mvmer`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md).

## Examples

``` r
if (.Platform$OS.type != "windows" || .Platform$r_arch !="i386") {
# \donttest{

#####
# Univariate joint model, with association structure based on the 
# current value of the linear predictor
f1 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year",
              # this next line is only to keep the example small in size!
              chains = 1, cores = 1, seed = 12345, iter = 1000)
print(f1) 
summary(f1) 
        
#####
# Univariate joint model, with association structure based on the 
# current value and slope of the linear predictor
f2 <- stan_jm(formulaLong = logBili ~ year + (year | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              assoc = c("etavalue", "etaslope"),
              time_var = "year",
              chains = 1, cores = 1, seed = 12345, iter = 1000)
print(f2)  

#####
# Univariate joint model, with association structure based on the 
# lagged value of the linear predictor, where the lag is 2 time 
# units (i.e. 2 years in this example)
f3 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year",
              assoc = "etavalue", lag_assoc = 2,
              chains = 1, cores = 1, seed = 12345, iter = 1000)
print(f3) 

#####
# Univariate joint model, where the association structure includes 
# interactions with observed data. Here we specify that we want to use 
# an association structure based on the current value of the linear 
# predictor from the longitudinal submodel (i.e. "etavalue"), but we 
# also want to interact this with the treatment covariate (trt) from
# pbcLong data frame, so that we can estimate a different association 
# parameter (i.e. estimated effect of log serum bilirubin on the log 
# hazard of death) for each treatment group
f4 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
              dataLong = pbcLong,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv,
              time_var = "year",
              assoc = c("etavalue", "etavalue_data(~ trt)"),
              chains = 1, cores = 1, seed = 12345, iter = 1000)
print(f4)

######
# Multivariate joint model, with association structure based 
# on the current value and slope of the linear predictor in the 
# first longitudinal submodel and the area under the marker 
# trajectory for the second longitudinal submodel
mv1 <- stan_jm(
        formulaLong = list(
          logBili ~ year + (1 | id), 
          albumin ~ sex + year + (year | id)),
        dataLong = pbcLong,
        formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
        dataEvent = pbcSurv,
        assoc = list(c("etavalue", "etaslope"), "etaauc"), 
        time_var = "year",
        chains = 1, cores = 1, seed = 12345, iter = 100)
print(mv1)

#####
# Multivariate joint model, where the association structure is formed by 
# including the expected value of each longitudinal marker (logBili and 
# albumin) in the linear predictor of the event submodel, as well as their 
# interaction effect (i.e. the interaction between the two "etavalue" terms). 
# Note that whether such an association structure based on a marker by 
# marker interaction term makes sense will depend on the context of your 
# application -- here we just show it for demostration purposes).
mv2 <- stan_jm(
        formulaLong = list(
          logBili ~ year + (1 | id), 
          albumin ~ sex + year + (year | id)),
        dataLong = pbcLong,
        formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
        dataEvent = pbcSurv,
        assoc = list(c("etavalue", "etavalue_etavalue(2)"), "etavalue"),
        time_var = "year", 
        chains = 1, cores = 1, seed = 12345, iter = 100)
        
#####
# Multivariate joint model, with one bernoulli marker and one
# Gaussian marker. We will artificially create the bernoulli
# marker by dichotomising log serum bilirubin
pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
mv3 <- stan_jm(
        formulaLong = list(
          ybern ~ year + (1 | id), 
          albumin ~ sex + year + (year | id)),
        dataLong = pbcLong,
        formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
        dataEvent = pbcSurv,
        family = list(binomial, gaussian),
        time_var = "year", 
        chains = 1, cores = 1, seed = 12345, iter = 1000)
# }
}
#> Fitting a univariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'jm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000197 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.97 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 2.475 seconds (Warm-up)
#> Chain 1:                1.689 seconds (Sampling)
#> Chain 1:                4.164 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> stan_jm
#>  formula (Long1): logBili ~ year + (1 | id)
#>  family  (Long1): gaussian [identity]
#>  formula (Event): Surv(futimeYears, death) ~ sex + trt
#>  baseline hazard: bs
#>  assoc:           etavalue (Long1)
#> ------
#> 
#> Longitudinal submodel: logBili
#>             Median MAD_SD
#> (Intercept) 0.802  0.203 
#> year        0.092  0.010 
#> sigma       0.509  0.026 
#> 
#> Event submodel:
#>                 Median MAD_SD exp(Median)
#> (Intercept)     -3.068  0.569  0.046     
#> sexf            -0.376  0.538  0.687     
#> trt             -0.719  0.506  0.487     
#> Long1|etavalue   1.449  0.281  4.257     
#> b-splines-coef1 -1.418  1.106     NA     
#> b-splines-coef2  0.222  0.859     NA     
#> b-splines-coef3 -1.768  1.192     NA     
#> b-splines-coef4  1.058  1.450     NA     
#> b-splines-coef5 -0.460  1.602     NA     
#> b-splines-coef6 -0.369  1.558     NA     
#> 
#> Group-level error terms:
#>  Groups Name              Std.Dev.
#>  id     Long1|(Intercept) 1.323   
#> Num. levels: id 40 
#> 
#> Sample avg. posterior predictive distribution 
#> of longitudinal outcomes:
#>                Median MAD_SD
#> Long1|mean_PPD 0.586  0.043 
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').Fitting a univariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'jm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000318 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.18 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 15.999 seconds (Warm-up)
#> Chain 1:                11.427 seconds (Sampling)
#> Chain 1:                27.426 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> stan_jm
#>  formula (Long1): logBili ~ year + (year | id)
#>  family  (Long1): gaussian [identity]
#>  formula (Event): Surv(futimeYears, death) ~ sex + trt
#>  baseline hazard: bs
#>  assoc:           etavalue (Long1), etaslope (Long1)
#> ------
#> 
#> Longitudinal submodel: logBili
#>             Median MAD_SD
#> (Intercept) 0.625  0.229 
#> year        0.249  0.063 
#> sigma       0.357  0.017 
#> 
#> Event submodel:
#>                 Median    MAD_SD    exp(Median)
#> (Intercept)        -3.416     0.757     0.033  
#> sexf               -0.440     0.720     0.644  
#> trt                -0.987     0.640     0.373  
#> Long1|etavalue      0.821     0.427     2.273  
#> Long1|etaslope      9.538     5.737 13874.955  
#> b-splines-coef1    -5.377     4.425        NA  
#> b-splines-coef2    -2.087     2.338        NA  
#> b-splines-coef3    -2.848     1.555        NA  
#> b-splines-coef4    -0.511     1.833        NA  
#> b-splines-coef5     0.326     1.747        NA  
#> b-splines-coef6    -0.537     1.694        NA  
#> 
#> Group-level error terms:
#>  Groups Name              Std.Dev. Corr
#>  id     Long1|(Intercept) 1.284        
#>         Long1|year        0.241    0.67
#> Num. levels: id 40 
#> 
#> Sample avg. posterior predictive distribution 
#> of longitudinal outcomes:
#>                Median MAD_SD
#> Long1|mean_PPD 0.585  0.027 
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').Fitting a univariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'jm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.00016 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.6 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 2.054 seconds (Warm-up)
#> Chain 1:                1.651 seconds (Sampling)
#> Chain 1:                3.705 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> stan_jm
#>  formula (Long1): logBili ~ year + (1 | id)
#>  family  (Long1): gaussian [identity]
#>  formula (Event): Surv(futimeYears, death) ~ sex + trt
#>  baseline hazard: bs
#>  assoc:           etavalue (Long1)
#> ------
#> 
#> Longitudinal submodel: logBili
#>             Median MAD_SD
#> (Intercept) 0.838  0.188 
#> year        0.092  0.010 
#> sigma       0.509  0.023 
#> 
#> Event submodel:
#>                 Median MAD_SD exp(Median)
#> (Intercept)     -2.910  0.631  0.054     
#> sexf            -0.380  0.591  0.684     
#> trt             -0.703  0.489  0.495     
#> Long1|etavalue   1.454  0.254  4.281     
#> b-splines-coef1 -1.600  1.023     NA     
#> b-splines-coef2  0.250  0.898     NA     
#> b-splines-coef3 -1.518  1.312     NA     
#> b-splines-coef4  0.913  1.666     NA     
#> b-splines-coef5 -0.029  1.747     NA     
#> b-splines-coef6 -0.174  1.602     NA     
#> 
#> Group-level error terms:
#>  Groups Name              Std.Dev.
#>  id     Long1|(Intercept) 1.325   
#> Num. levels: id 40 
#> 
#> Sample avg. posterior predictive distribution 
#> of longitudinal outcomes:
#>                Median MAD_SD
#> Long1|mean_PPD 0.585  0.042 
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').Fitting a univariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'jm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000174 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.74 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 3.098 seconds (Warm-up)
#> Chain 1:                2.201 seconds (Sampling)
#> Chain 1:                5.299 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.06, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> stan_jm
#>  formula (Long1): logBili ~ year + (1 | id)
#>  family  (Long1): gaussian [identity]
#>  formula (Event): Surv(futimeYears, death) ~ sex + trt
#>  baseline hazard: bs
#>  assoc:           etavalue (Long1), etavalue_data (Long1)
#> ------
#> 
#> Longitudinal submodel: logBili
#>             Median MAD_SD
#> (Intercept) 0.830  0.206 
#> year        0.092  0.008 
#> sigma       0.508  0.022 
#> 
#> Event submodel:
#>                    Median MAD_SD exp(Median)
#> (Intercept)        -3.140  0.690  0.043     
#> sexf               -0.366  0.629  0.694     
#> trt                -0.184  0.929  0.832     
#> Long1|etavalue      1.598  0.323  4.945     
#> Long1|etavalue:trt -0.374  0.574  0.688     
#> b-splines-coef1    -1.477  1.254     NA     
#> b-splines-coef2     0.038  0.856     NA     
#> b-splines-coef3    -1.639  1.234     NA     
#> b-splines-coef4     0.736  1.674     NA     
#> b-splines-coef5    -0.292  1.705     NA     
#> b-splines-coef6    -0.227  1.665     NA     
#> 
#> Group-level error terms:
#>  Groups Name              Std.Dev.
#>  id     Long1|(Intercept) 1.345   
#> Num. levels: id 40 
#> 
#> Sample avg. posterior predictive distribution 
#> of longitudinal outcomes:
#>                Median MAD_SD
#> Long1|mean_PPD 0.588  0.039 
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').Fitting a multivariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'jm' NOW (CHAIN 1).
#> Chain 1: Rejecting initial value:
#> Chain 1:   Log probability evaluates to log(0), i.e. negative infinity.
#> Chain 1:   Stan can't start sampling from this initial value.
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.001553 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 15.53 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 7
#> Chain 1:            adapt_window = 38
#> Chain 1:            term_buffer = 5
#> Chain 1: 
#> Chain 1: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 1: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 1: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 1: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 1: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 1: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 1: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 1: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 1: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 1: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 1: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 1: Iteration: 100 / 100 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 18.655 seconds (Warm-up)
#> Chain 1:                0.509 seconds (Sampling)
#> Chain 1:                19.164 seconds (Total)
#> Chain 1: 
#> Warning: There were 50 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 2.15, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Markov chains did not converge! Do not analyze results!
#> stan_jm
#>  formula (Long1): logBili ~ year + (1 | id)
#>  family  (Long1): gaussian [identity]
#>  formula (Long2): albumin ~ sex + year + (year | id)
#>  family  (Long2): gaussian [identity]
#>  formula (Event): Surv(futimeYears, death) ~ sex + trt
#>  baseline hazard: bs
#>  assoc:           etavalue (Long1), etaslope (Long1), etaauc (Long2)
#> ------
#> 
#> Longitudinal submodel 1: logBili
#>             Median MAD_SD
#> (Intercept) 0.697  0.000 
#> year        0.086  0.000 
#> sigma       0.526  0.000 
#> 
#> Longitudinal submodel 2: albumin
#>             Median MAD_SD
#> (Intercept)  3.509  0.000
#> sexf         0.081  0.000
#> year        -0.123  0.000
#> sigma        0.334  0.000
#> 
#> Event submodel:
#>                 Median        MAD_SD        exp(Median)  
#> (Intercept)      1.082885e+11  1.330000e-01           Inf
#> sexf            -3.470000e-01  0.000000e+00  7.070000e-01
#> trt             -9.600000e-02  0.000000e+00  9.080000e-01
#> Long1|etavalue   3.422000e+00  0.000000e+00  3.064300e+01
#> Long1|etaslope  -1.264015e+12  1.549000e+00  0.000000e+00
#> Long2|etaauc     3.360000e-01  0.000000e+00  1.399000e+00
#> b-splines-coef1  0.000000e+00  0.000000e+00            NA
#> b-splines-coef2  0.000000e+00  0.000000e+00            NA
#> b-splines-coef3  0.000000e+00  0.000000e+00            NA
#> b-splines-coef4  0.000000e+00  0.000000e+00            NA
#> b-splines-coef5  0.000000e+00  0.000000e+00            NA
#> b-splines-coef6  0.000000e+00  0.000000e+00            NA
#> 
#> Group-level error terms:
#>  Groups Name              Std.Dev. Corr     
#>  id     Long1|(Intercept) 1.11688           
#>         Long2|(Intercept) 0.36463  0.00     
#>         Long2|year        0.03419  0.00 0.00
#> Num. levels: id 40 
#> 
#> Sample avg. posterior predictive distribution 
#> of longitudinal outcomes:
#>                Median MAD_SD
#> Long1|mean_PPD 0.885  0.036 
#> Long2|mean_PPD 3.145  0.021 
#> 
#> ------
#> For info on the priors used see help('prior_summary.stanreg').Fitting a multivariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'jm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000333 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.33 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 7
#> Chain 1:            adapt_window = 38
#> Chain 1:            term_buffer = 5
#> Chain 1: 
#> Chain 1: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 1: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 1: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 1: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 1: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 1: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 1: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 1: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 1: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 1: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 1: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 1: Iteration: 100 / 100 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 2.73 seconds (Warm-up)
#> Chain 1:                6.373 seconds (Sampling)
#> Chain 1:                9.103 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.14, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Fitting a multivariate joint model.
#> 
#> Please note the warmup may be much slower than later iterations!
#> 
#> SAMPLING FOR MODEL 'jm' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000257 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.57 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 12.701 seconds (Warm-up)
#> Chain 1:                5.811 seconds (Sampling)
#> Chain 1:                18.512 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.06, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
```
