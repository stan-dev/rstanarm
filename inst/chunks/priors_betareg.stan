  // Log-priors for coefficients
  if (prior_dist_z == 1)  target += normal_lpdf(z_omega | 0, 1);
  else if (prior_dist_z == 2) target += normal_lpdf(z_omega | 0, 1);
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
  else if (prior_dist_z == 5) { // laplace
    target += normal_lpdf(z_omega | 0, 1);
    target += exponential_lpdf(S_z[1] | 1);
  }
  else if (prior_dist_z == 6) { // lasso
    target += normal_lpdf(z_omega | 0, 1);
    target += exponential_lpdf(S_z[1] | 1);
    target += chi_square_lpdf(one_over_lambda_z[1] | prior_df_z[1]);
  }
  else if (prior_dist_z == 7) { // product_normal
    target += normal_lpdf(z_omega | 0, 1);
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
