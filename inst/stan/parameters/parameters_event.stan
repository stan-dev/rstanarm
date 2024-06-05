  array[e_has_intercept] real e_gamma; // intercept for event submodel
  vector[e_K] e_z_beta;          // primitive log hazard ratios

  // unscaled basehaz params, either:
  //   - weibull shape parameter
  //   - b-spline coefs on log basehaz
  //   - coefs for piecewise constant basehaz
  vector<lower=(basehaz_type == 1 ? 0 : negative_infinity())>[basehaz_df] e_aux_unscaled;

  // parameters for priors on log haz ratios
  array[e_hs] real<lower=0> e_global;
  array[e_hs] vector<lower=0>[(e_hs>0)*e_K] e_local;
  array[e_hs > 0] real<lower=0> e_caux;
  array[e_prior_dist == 5 || e_prior_dist == 6] vector<lower=0>[e_K] e_mix;
  array[e_prior_dist == 6] real<lower=0> e_ool;
