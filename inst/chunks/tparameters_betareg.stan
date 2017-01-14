  if (prior_dist_z == 0) omega = z_omega;
  else if (prior_dist_z == 1) omega = z_omega .* prior_scale_z + prior_mean_z;
  else if (prior_dist_z == 2) for (k in 1:z_dim) {
    omega[k] = CFt(omega[k], prior_df_z[k]) * prior_scale_z[k] + prior_mean_z[k];
  }
  else if (prior_dist_z == 3) omega = hs_prior(z_omega, global_z, local_z);
  else if (prior_dist_z == 4) omega = hsplus_prior(z_omega, global_z, local_z);
