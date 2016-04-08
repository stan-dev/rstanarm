  // Log-priors for coefficients
  if (prior_dist == 1) z_beta ~ normal(0, 1);
  else if (prior_dist == 2) {
    if (t_all_124) z_beta ~ normal(0,1);
    else if (t_any_124) for (k in 1:K) {
      if (prior_df[k] == 1 || prior_df[k] == 2 || prior_df[k] == 4)
        z_beta[k] ~ normal(0,1);
      else z_beta[k] ~ student_t(prior_df[k], 0, 1);
    }
    else z_beta ~ student_t(prior_df, 0, 1);
  }
  else if (prior_dist == 3) { // hs
    z_beta ~ normal(0,1);
    local[1] ~ normal(0,1);
    local[2] ~ inv_gamma(0.5 * prior_df, 0.5 * prior_df);
    global[1] ~ normal(0,1);
    global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (prior_dist == 4) { // hs+
    z_beta ~ normal(0,1);
    local[1] ~ normal(0,1);
    local[2] ~ inv_gamma(0.5 * prior_df, 0.5 * prior_df);
    local[3] ~ normal(0,1);
    // unorthodox useage of prior_scale as another df hyperparameter
    local[4] ~ inv_gamma(0.5 * prior_scale, 0.5 * prior_scale);
    global[1] ~ normal(0,1);
    global[2] ~ inv_gamma(0.5, 0.5);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1)  // normal
      gamma ~ normal(prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2)  // student_t
      gamma ~ student_t(prior_df_for_intercept, prior_mean_for_intercept, 
                        prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }
