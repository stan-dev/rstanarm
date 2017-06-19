  vector[prior_dist == 7 ? sum(num_normals) : K] z_beta;
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  vector<lower=0>[K] S[prior_dist == 5 || prior_dist == 6];
  real<lower=0> one_over_lambda[prior_dist == 6];
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t - (special_case && interaction_prior == 1 ? n_multi_way : 0)] tau;
  vector<lower=0>[len_multi_way_uniq] lambda_multi_way[special_case == 1
                                                       && n_multi_way > 0
                                                       && interaction_prior == 1];
  real<lower=0> glob_scale[n_multi_way > 0 && interaction_prior == 1 && special_case == 1];
