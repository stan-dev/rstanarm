  // hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  vector<lower=0>[M] prior_scale_for_intercept;
  vector<lower=0>[M] prior_scale_for_aux;
  vector[K]          prior_mean;
  vector[M]          prior_mean_for_intercept;
  vector<lower=0>[M] prior_mean_for_aux;
  vector<lower=0>[K] prior_df;
  vector<lower=0>[M] prior_df_for_intercept;
  vector<lower=0>[M] prior_df_for_aux;
  vector<lower=0>[M] global_prior_df;    // for hs priors only 
  vector<lower=0>[M] global_prior_scale; // for hs priors only
  
  // 1 = same prior type for all submodels (none, normal or student-t)
  int<lower=0,upper=1> prior_special_case; 
