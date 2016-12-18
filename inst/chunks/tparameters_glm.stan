  vector[K] beta;
  vector[q] b;
  vector[len_theta_L] theta_L;
  if      (prior_dist == 0) beta = z_beta;
  else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 3) 
    beta = hs_prior(z_beta, global, local, global_prior_scale);
  else if (prior_dist == 4) 
    beta = hsplus_prior(z_beta, global, local);
