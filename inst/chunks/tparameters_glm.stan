  vector[K] beta;
  vector[K_smooth] beta_smooth;
  vector[K_smooth > 0 ? smooth_map[K_smooth] : 0] smooth_sd;
  vector[q] b;
  vector[len_theta_L] theta_L;

  #include "tparameters.stan"

  if (K_smooth) {
    smooth_sd = prior_mean_for_smooth + prior_scale_for_smooth .* smooth_sd_raw;
    if (is_continuous && family == 1) smooth_sd = smooth_sd * aux;
    beta_smooth = z_beta_smooth .* smooth_sd[smooth_map];
  }
