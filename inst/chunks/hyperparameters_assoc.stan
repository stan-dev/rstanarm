  // hyperparameter values are set to 0 if there is no prior
  vector[a_K]          a_prior_mean;
  vector<lower=0>[a_K] a_prior_scale;
  vector<lower=0>[a_K] a_prior_df;
  real<lower=0>        a_global_prior_scale; // for hs priors only
  real<lower=0>        a_global_prior_df;
  real<lower=0>        a_slab_df;
  real<lower=0>        a_slab_scale;
