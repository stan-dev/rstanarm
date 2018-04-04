  // hyperparameter values are set to 0 if there is no prior
  vector[e_K]                 e_prior_mean;
  real                        e_prior_mean_for_intercept;
  vector[basehaz_df]          e_prior_mean_for_aux;
  vector<lower=0>[e_K]        e_prior_scale;
  real<lower=0>               e_prior_scale_for_intercept;
  vector<lower=0>[basehaz_df] e_prior_scale_for_aux;
  vector<lower=0>[e_K]        e_prior_df;
  real<lower=0>               e_prior_df_for_intercept;
  vector<lower=0>[basehaz_df] e_prior_df_for_aux;
  real<lower=0>               e_global_prior_scale; // for hs priors only
  real<lower=0>               e_global_prior_df;
  real<lower=0>               e_slab_df;
  real<lower=0>               e_slab_scale;
