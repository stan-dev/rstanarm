  real poisson_max = pow(2.0, 30.0);
  real sum_log_y[M];
  vector[N_real] sqrt_y = rep_vector(not_a_number(), N_real);
  vector[N_real] log_y  = rep_vector(not_a_number(), N_real);
  int<lower=0> hsM[M]; // used instead of hs... hs from tdata_glm.stan will be ignored
  int<lower=0> len_global = 0;  // total num. hs global params
  int<lower=0> len_local2 = 0;  // total num. hs local params
  int<lower=0> len_local4 = 0;  // total num. hs_plus local params
  int<lower=0> len_mix    = 0;  // total num. shrinkage params
  int<lower=0> len_ool    = 0;  // total num. one over lambda params
  int<lower=0> len_noise  = 0;  // total num. noise params
  int<lower=0> idx_global[M,2]; // indexing for hs global params
  int<lower=0> idx_local2[M,2]; // indexing for hs local params
  int<lower=0> idx_local4[M,2]; // indexing for hs_plus local params
  int<lower=0> idx_mix[M,2];    // indexing for shrinkage params
  int<lower=0> idx_ool[M];      // indexing for one over lambda params
  int<lower=0> idx_noise[M,2];  // indexing for noise params
  vector[0] beta_smooth;        // not used

  // the following is almost the same as tdata_glm.stan
  // but can't use that file because prior_dist throws error 
  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p) * (t > 0);
  int<lower=0> len_rho = sum(p) - t;
  int<lower=0, upper=1> is_continuous = 0; // changed in continuous.stan
  int<lower=1> pos = 1;
  real<lower=0> delta[len_concentration];  
