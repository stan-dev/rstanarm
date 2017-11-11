  vector[a_K] a_z_beta; // primitive assoc params

  // parameters for priors on assoc params 
  real<lower=0> a_global[a_hs];
  vector<lower=0>[(a_hs>0)*a_K] a_local[a_hs];
  vector<lower=0>[a_K] a_mix[a_prior_dist == 5 || a_prior_dist == 6];
  real<lower=0> a_ool[a_prior_dist == 6];
