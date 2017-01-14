  vector[z_dim] z_omega;          // betareg z variable coefficients
  real gamma_z[has_intercept_z];  // betareg intercept
  real<lower=0> global_z[hs_z];
  vector<lower=0>[z_dim] local_z[hs_z];
