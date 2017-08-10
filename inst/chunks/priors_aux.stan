// Log-priors
if (prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
  real log_half = -0.693147180559945286;
  if (prior_dist_for_aux == 1)
    target += normal_lpdf(aux_unscaled | 0, 1) - log_half;
  else if (prior_dist_for_aux == 2)
    target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1) - log_half;
  else
   target += exponential_lpdf(aux_unscaled | 1);
}
