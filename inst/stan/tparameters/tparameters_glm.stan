  vector[K] beta;
  vector[K_smooth] beta_smooth;
  vector[K_smooth > 0 ? smooth_map[K_smooth] : 0] smooth_sd;
  vector[q] b;
  vector[len_theta_L] theta_L;
  if      (prior_dist == 0) beta = z_beta;
  else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 3) {
    real c2 = square(slab_scale) * caux[1];
    if (is_continuous == 1 && family == 1)
      beta = hs_prior(z_beta, global, local, global_prior_scale, aux, c2);
    else beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 4) {
    real c2 = square(slab_scale) * caux[1];
    if (is_continuous == 1 && family == 1)
      beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux, c2);
    else beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 5) // laplace
    beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 6) // lasso
    beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 7) { // product_normal
    int z_pos = 1;
    for (k in 1:K) {
      beta[k] = z_beta[z_pos];
      z_pos += 1;
      for (n in 2:num_normals[k]) {
        beta[k] *= z_beta[z_pos];
        z_pos += 1;
      }
      beta[k] *= prior_scale[k] ^ num_normals[k];
      beta[k] += prior_mean[k];
    }
  }

  if (K_smooth) {
    smooth_sd = prior_mean_for_smooth + prior_scale_for_smooth .* smooth_sd_raw;
    if (is_continuous && family == 1) smooth_sd *= aux;
    beta_smooth = z_beta_smooth .* smooth_sd[smooth_map];
  }
