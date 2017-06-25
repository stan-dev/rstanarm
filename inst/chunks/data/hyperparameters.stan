  // hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_scale_for_aux;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  real<lower=0> prior_mean_for_aux;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_df_for_aux;  
  real<lower=0> global_prior_df;    // for hs priors only
  real<lower=0> global_prior_scale; // for hs priors only
  int<lower=2> num_normals[prior_dist == 7 ? K : 0];
