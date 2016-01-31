  // hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_scale_for_dispersion;
