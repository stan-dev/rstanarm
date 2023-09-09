  vector[prior_dist_z == 7 ? sum(num_normals_z) : z_dim] z_omega; // betareg z variable coefficients
  array[has_intercept_z] real<lower=(link_phi <= 1 ? negative_infinity() : 0)> gamma_z;  // betareg intercept
  array[hs_z] real<lower=0> global_z;
  array[hs_z] vector<lower=0>[z_dim] local_z;
  array[hs_z > 0] real<lower=0> caux_z;
  array[prior_dist_z == 5 || prior_dist_z == 6] vector<lower=0>[z_dim] S_z;
  array[prior_dist_z == 6] real<lower=0> one_over_lambda_z;
