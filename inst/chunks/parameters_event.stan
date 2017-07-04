  real e_gamma[e_has_intercept]; // intercept for event submodel 
  vector[e_K] e_z_beta;          // primitive log hazard ratios
  
  // unscaled basehaz params, either:
  //   - weibull shape parameter
  //   - b-spline coefs on log basehaz
  //   - coefs for piecewise constant basehaz
  vector<lower=(basehaz_type == 1 ? 0 : negative_infinity())>[basehaz_df] e_aux_unscaled;       

  // parameters for priors on log haz ratios
  real<lower=0> e_global[e_hs];
  vector<lower=0>[(e_hs>0)*e_K] e_local[e_hs];
  vector<lower=0>[e_K] e_mix[e_prior_dist == 5 || e_prior_dist == 6];
  real<lower=0> e_ool[e_prior_dist == 6];
