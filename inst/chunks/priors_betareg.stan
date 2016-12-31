  // Log-priors for coefficients
  if (prior_dist_z == 1)  target += normal_lpdf(z_omega | 0, 1);
  else if (prior_dist_z == 2) {
    if (t_all_124_z) target += normal_lpdf(z_omega | 0, 1);
    else if (t_any_124_z) for (k in 1:z_dim) {
      if (prior_df_z[k] == 1 || prior_df_z[k] == 2 || prior_df_z[k] == 4)
        target += normal_lpdf(z_omega[k] | 0,1);
      else target += student_t_lpdf(z_omega[k] | prior_df_z[k], 0, 1);
    }
    else target += student_t_lpdf(z_omega | prior_df_z, 0, 1);
  }
  else if (prior_dist_z == 3) { // hs
    target += normal_lpdf(z_omega | 0, 1);
    target += normal_lpdf(local_z[1] | 0, 1);
    target += inv_gamma_lpdf(local_z[2] | 0.5 * prior_df_z, 0.5 * prior_df_z);
    target += normal_lpdf(global_z[1] | 0, 1);
    target += inv_gamma_lpdf(global_z[2] | 0.5, 0.5);
  }
  else if (prior_dist_z == 4) { // hs+
    target += normal_lpdf(z_omega | 0, 1);
    target += normal_lpdf(local_z[1] | 0, 1);
    target += inv_gamma_lpdf(local_z[2] | 0.5 * prior_df_z, 0.5 * prior_df_z);
    target += normal_lpdf(local_z[3] | 0, 1);
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local_z[4] | 0.5 * prior_scale_z, 0.5 * prior_scale_z);
    target += normal_lpdf(global_z[1] | 0, 1);
    target += inv_gamma_lpdf(global_z[2] | 0.5, 0.5);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept_z == 1) {
    if (prior_dist_for_intercept_z == 1)  // normal
      target += normal_lpdf(gamma_z | prior_mean_for_intercept_z, prior_scale_for_intercept_z);
    else if (prior_dist_for_intercept_z == 2)  // student_t
      target += student_t_lpdf(gamma_z | prior_df_for_intercept_z, prior_mean_for_intercept_z, 
                               prior_scale_for_intercept_z);
    /* else prior_dist is 0 and nothing is added */
  }
