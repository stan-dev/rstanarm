  vector[K] beta;
  vector[q] b;
  vector[len_theta_L] theta_L;
  if      (prior_dist == 0) beta = z_beta;
  else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 3) beta = hs_prior(z_beta, global, local);
  else if (prior_dist == 4) beta = hsplus_prior(z_beta, global, local);
  else if (prior_dist == 5) // laplace
    beta = prior_mean + prior_scale .* sqrt(2 * V[1]) .* z_beta;
  else if (prior_dist == 6) // lasso
    beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * V[1]) .* z_beta;
  else if (prior_dist == 7) { // product_normal
    int z_pos = 1;
    for (k in 1:K) {
      vector[num_normals[k]] ord;
      ord[1] = z_beta[z_pos];
      z_pos = z_pos + 1;
      beta[k] = 1;
      for (n in 2:rows(ord)) {
        real z = z_beta[z_pos];
        ord[n] = ord[n - 1] + exp(z); // Jacobian in priors_glm.stan
        beta[k] = beta[k] * ord[n];
        z_pos = z_pos + 1;
      }
      beta[k] = beta[k] * prior_scale[k] ^ rows(ord) + prior_mean[k];
    }
    if (pos != rows(z_beta) + 1) reject("wrong number of parameters");
  }
