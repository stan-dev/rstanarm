# Fitted model objects

The rstanarm model-fitting functions return an object of class
`'stanreg'`, which is a list containing at a minimum the components
listed below. Each `stanreg` object will also have additional classes
(e.g. 'aov', 'betareg', 'glm', 'polr', etc.) and several additional
components depending on the model and estimation algorithm.  
  
Some additional details apply to models estimated using the
[`stan_mvmer`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md) or
[`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md) modelling
functions. The
[`stan_mvmer`](https://mc-stan.org/rstanarm/reference/stan_mvmer.md)
modelling function returns an object of class `'stanmvreg'`, which
inherits the `'stanreg'` class, but has a number of additional elements
described in the subsection below. The
[`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md) modelling
function returns an object of class `'stanjm'`, which inherits both the
`'stanmvreg'` and `'stanreg'` classes, but has a number of additional
elements described in the subsection below. Both the `'stanjm'` and
`'stanmvreg'` classes have several of their own methods for situations
in which the default `'stanreg'` methods are not suitable; see the **See
Also** section below.

## Note

The [`stan_biglm`](https://mc-stan.org/rstanarm/reference/stan_biglm.md)
function is an exception. It returns a
[stanfit](https://mc-stan.org/rstan/reference/stanfit-class.html) object
rather than a stanreg object.

## Elements for `stanreg` objects

- `coefficients`:

  Point estimates, as described in
  [`print.stanreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md).

- `ses`:

  Standard errors based on [`mad`](https://rdrr.io/r/stats/mad.html), as
  described in
  [`print.stanreg`](https://mc-stan.org/rstanarm/reference/print.stanreg.md).

- `residuals`:

  Residuals of type `'response'`.

- `fitted.values`:

  Fitted mean values. For GLMs the linear predictors are transformed by
  the inverse link function.

- `linear.predictors`:

  Linear fit on the link scale. For linear models this is the same as
  `fitted.values`.

- `covmat`:

  Variance-covariance matrix for the coefficients based on draws from
  the posterior distribution, the variational approximation, or the
  asymptotic sampling distribution, depending on the estimation
  algorithm.

- `model,x,y`:

  If requested, the the model frame, model matrix and response variable
  used, respectively.

- `family`:

  The [`family`](https://rdrr.io/r/stats/family.html) object used.

- `call`:

  The matched call.

- `formula`:

  The model [`formula`](https://rdrr.io/r/stats/formula.html).

- `data,offset,weights`:

  The `data`, `offset`, and `weights` arguments.

- `algorithm`:

  The estimation method used.

- `prior.info`:

  A list with information about the prior distributions used.

- `stanfit,stan_summary`:

  The object of
  [`stanfit-class`](https://mc-stan.org/rstan/reference/stanfit-class.html)
  returned by RStan and a matrix of various summary statistics from the
  stanfit object.

- `rstan_version`:

  The version of the rstan package that was used to fit the model.

## Elements for `stanmvreg` objects

`cnms`

:   The names of the grouping factors and group specific parameters,
    collapsed across the longitudinal or glmer submodels.

`flevels`

:   The unique factor levels for each grouping factor, collapsed across
    the longitudinal or glmer submodels.

`n_markers`

:   The number of longitudinal or glmer submodels.

`n_yobs`

:   The number of observations for each longitudinal or glmer submodel.

`n_grps`

:   The number of levels for each grouping factor (for models estimated
    using
    [`stan_jm`](https://mc-stan.org/rstanarm/reference/stan_jm.md), this
    will be equal to `n_subjects` if the individual is the only grouping
    factor).

`runtime`

:   The time taken to fit the model (in minutes).

## Additional elements for `stanjm` objects

- `id_var,time_var`:

  The names of the variables distinguishing between individuals, and
  representing time in the longitudinal submodel.

- `n_subjects`:

  The number of individuals.

- `n_events`:

  The number of non-censored events.

- `eventtime,status`:

  The event (or censoring) time and status indicator for each
  individual.

- `basehaz`:

  A list containing information about the baseline hazard.

- `assoc`:

  An array containing information about the association structure.

- `epsilon`:

  The width of the one-sided difference used to numerically evaluate the
  slope of the longitudinal trajectory; only relevant if a slope-based
  association structure was specified (e.g. etaslope, muslope, etc).

- `qnodes`:

  The number of Gauss-Kronrod quadrature nodes used to evaluate the
  cumulative hazard in the joint likelihood function.

## See also

[`stanreg-methods`](https://mc-stan.org/rstanarm/reference/stanreg-methods.md),
[`stanmvreg-methods`](https://mc-stan.org/rstanarm/reference/stanmvreg-methods.md)
