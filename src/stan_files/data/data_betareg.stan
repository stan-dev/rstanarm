  // betareg data
  int<lower=0, upper=1> has_intercept_z;  // presence of z intercept
  int<lower=0> link_phi;                  // link transformation for eta_z (0 => no z in model)
  int<lower=0> z_dim;                     // dimensions of z vars
  matrix[N, z_dim] betareg_z;             // matrix of z vars
  row_vector[z_dim] zbar;                 // mean of predictors
  // betareg hyperparameters
  int<lower=0,upper=7> prior_dist_z;
  int<lower=0,upper=2> prior_dist_for_intercept_z;
  vector<lower=0>[z_dim] prior_scale_z;
  real<lower=0> prior_scale_for_intercept_z;
  vector[z_dim] prior_mean_z;
  real prior_mean_for_intercept_z;
  vector<lower=0>[z_dim] prior_df_z;
  real<lower=0> prior_df_for_intercept_z;
  real<lower=0> global_prior_scale_z;
  real<lower=0> global_prior_df_z;
  real<lower=0> slab_df_z;
  real<lower=0> slab_scale_z;
  array[prior_dist_z == 7 ? z_dim : 0] int<lower=2> num_normals_z;
