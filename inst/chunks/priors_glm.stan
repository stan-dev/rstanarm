  #include "priors.stan"
  
  if (K_smooth) {
    target += normal_lpdf(z_beta_smooth | 0, 1);
    if (prior_dist_for_smooth > 0) {
      real log_half = -0.693147180559945286;
      if (prior_dist_for_smooth == 1)
        target += normal_lpdf(smooth_sd_raw | 0, 1) - log_half;
      else if (prior_dist_for_smooth == 2)
        target += student_t_lpdf(smooth_sd_raw | prior_df_for_smooth, 0, 1) - log_half;
      else if (prior_dist_for_smooth == 3)
        target += exponential_lpdf(smooth_sd_raw | 1);
    }
  }
