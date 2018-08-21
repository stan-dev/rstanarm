  if (prior_dist_z == 0) omega = z_omega;
  else if (prior_dist_z == 1) omega = z_omega .* prior_scale_z + prior_mean_z;
  else if (prior_dist_z == 2) for (k in 1:z_dim) {
    real left = CFt(omega[k], prior_df_z[k]); 
    omega[k] = left * prior_scale_z[k] + prior_mean_z[k];
  }
  else if (prior_dist_z == 3) 
    omega = hs_prior(z_omega, global_z, local_z, global_prior_scale, 
                     1, square(slab_scale_z) * caux_z[1]);
  else if (prior_dist_z == 4) 
    omega = hsplus_prior(z_omega, global_z, local_z, global_prior_scale, 1,
                         square(slab_scale_z) * caux_z[1]);
  else if (prior_dist_z == 5)
    omega = prior_mean_z + prior_scale_z .* sqrt(2 * S_z[1]) .* z_omega;
  else if (prior_dist_z == 6)
    omega = prior_mean_z + one_over_lambda_z[1] * prior_scale_z .* sqrt(2 * S_z[1]) .* z_omega;
  else if (prior_dist_z == 7) {
    int z_pos = 1;
    for (k in 1:z_dim) {
      omega[k] = z_omega[z_pos];
      z_pos += 1;
      for (n in 2:num_normals_z[k]) {
        omega[k] *= z_omega[z_pos];
        z_pos += 1;
      }
      omega[k] *= prior_scale_z[k] ^ num_normals_z[k];
      omega[k] += prior_mean_z[k];
    }
  }
    
