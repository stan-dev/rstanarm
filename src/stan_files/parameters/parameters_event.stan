  // primitive log hazard ratios
  vector[e_K] e_z_beta;

  // intercept
  real e_gamma[e_has_intercept == 1];

  // unscaled basehaz parameters
  //   exp model:      nvars = 0, ie. no aux parameter
  //   weibull model:  nvars = 1, ie. shape parameter
  //   gompertz model: nvars = 1, ie. scale parameter
  //   M-spline model: nvars = number of basis terms, ie. spline coefs
  //   B-spline model: nvars = number of basis terms, ie. spline coefs
  vector<lower=coefs_lb(basehaz_type)>[basehaz_nvars] e_aux_unscaled;

  // parameters for priors on log haz ratios
  real<lower=0> e_global[e_hs];
  vector<lower=0>[(e_hs>0)*e_K] e_local[e_hs];
  real<lower=0> e_caux[e_hs > 0];
  vector<lower=0>[e_K] e_mix[e_prior_dist == 5 || e_prior_dist == 6];
  real<lower=0> e_ool[e_prior_dist == 6];
