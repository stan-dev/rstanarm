  // the last section of this file is the same as parameters_glm.stan,
  // therefore that file could be split in two to save duplication
  // (e.g. split into "parameters_priors.stan" and "parameters_glmer.stan")

  // intercepts, differs from parameters_glm.stan
  real          gamma_nob[sum_has_intercept_nob]; // unbounded
  real<lower=0> gamma_lob[sum_has_intercept_lob]; // lower bounded 
  real<upper=0> gamma_upb[sum_has_intercept_upb]; // upper bounded 
  
  // differs from first part of parameters_glm.stan
  vector[K] z_beta;
  real<lower=0> global[len_global];
  vector<lower=0>[len_local2] local2[(len_local2 > 0) ? 2 : 0];
  vector<lower=0>[len_local4] local4[(len_local4 > 0) ? 4 : 0];
  vector<lower=0>[len_mix] mix[(len_mix > 0)];
  real<lower=0> ool[len_ool]; // one_over_lambda
  vector<lower=0>[len_noise] noise[(len_noise > 0)];
  vector<lower=0>[sum_has_aux] aux_unscaled; // interpretation depends on family!
  
  // same as second part of parameters_glm.stan
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
