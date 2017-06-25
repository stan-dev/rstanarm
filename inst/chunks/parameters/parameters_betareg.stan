  vector[prior_dist_z == 7 ? sum(num_normals_z) : z_dim] z_omega; // betareg z variable coefficients
  real<lower=(link_phi <= 1 ? negative_infinity() : 0)> gamma_z[has_intercept_z];  // betareg intercept
  real<lower=0> global_z[hs_z];
  vector<lower=0>[z_dim] local_z[hs_z];
  vector<lower=0>[z_dim] S_z[prior_dist_z == 5 || prior_dist_z == 6];
  real<lower=0> one_over_lambda_z[prior_dist_z == 6];
