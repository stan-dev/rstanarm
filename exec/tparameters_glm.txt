  vector[K] beta;
  vector[q] b;
  vector[len_theta_L] theta_L;
  if      (prior_dist == 0) beta <- z_beta;
  else if (prior_dist <= 2) beta <- z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 3) beta <- hs_prior(z_beta, global, local);
  else if (prior_dist == 4) beta <- hsplus_prior(z_beta, global, local);
