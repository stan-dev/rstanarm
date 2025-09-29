  vector[a_K] a_z_beta; // primitive assoc params

  // parameters for priors on assoc params
  array[a_hs] real<lower=0> a_global;
  array[a_hs] vector<lower=0>[(a_hs>0)*a_K] a_local;
  array[a_hs > 0] real<lower=0> a_caux;
  array[a_prior_dist == 5 || a_prior_dist == 6] vector<lower=0>[a_K] a_mix;
  array[a_prior_dist == 6] real<lower=0> a_ool;
