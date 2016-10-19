  if (prior_dist_z == 0) omega = z_omega;
  else if (prior_dist_z == 1) omega = z_omega .* prior_scale_z + prior_mean_z;
  else if (prior_dist_z == 2) for (k in 1:z_dim) {
    real P_z;
    if (prior_df_z[k] == 1) {
      P_z = Phi(z_omega[k]);
      omega[k] = tan(pi() * (P_z - 0.5));
    }
    else if (prior_df_z[k] == 2) {
      P_z = Phi(z_omega[k]);
      omega[k] = 2 * (P_z - 0.5) / sqrt(2.0 * P_z * (1 - P_z));
    }
    else if (prior_df_z[k] == 4) {
      real q_a_z;
      P_z = Phi(z_omega[k]);
      q_a_z = sqrt(4.0 * P_z * (1 - P_z));
      q_a_z = cos(acos(q_a_z) / 3) / q_a_z;
      omega[k] = 2 * sqrt(q_a_z - 1);
      if (P_z < 0.5) omega[k] = -omega[k];
    }
    else omega[k] = z_omega[k];
    omega[k] = omega[k] * prior_scale_z[k] + prior_mean_z[k];
  }
  else if (prior_dist_z == 3) omega = hs_prior(z_omega, global_z, local_z);
  else if (prior_dist_z == 4) omega = hsplus_prior(z_omega, global_z, local_z);
