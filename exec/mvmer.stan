#include "Columbia_copyright.stan"
#include "Brilleman_copyright.stan"
#include "license.stan" // GPL3+

// Multivariate GLM with group-specific terms
functions {
  #include "common_functions.stan"
  #include "bernoulli_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "continuous_likelihoods.stan"
  #include "count_likelihoods.stan"
  #include "jm_functions.stan"
}
data {
  #include "NKX.stan" // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x  
  
  // declares M, N{M,_real,_int,01}, idx{_real,_int,_K,_hs2,hs4}, KM,
  //   sum_has_intercept{_nob,_lob,_upb}, sum_has_aux
  #include "dimensions_mvglm.stan" 
  
  // declares prior_PD, has_intercept{_nob,_lob,_upb}, family, link, prior_dist, 
  //   prior_dist_for_intercept, prior_dist_for_aux
  #include "data_mvglm.stan"     // same as data_glm.stan, but arrays of size M 
  #include "data2_mvglm.stan"    // declares y_{real,int}, has_aux 
  #include "weights_offset.stan" // declares has_weights, weights, has_offset, offset
 
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, 
  //   {len_}concentration, {len_}regularization
  #include "glmer_stuff.stan"  
  #include "glmer_stuff2.stan" // declares num_not_zero, w, v, u, special_case
  #include "mvmer_stuff.stan"  // declares pmat, qmat, q1, q2

  // declares {e_,a_}{prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, 
  //   prior_{mean, scale, df}_for_aux, global_prior_{df,scale}}
  #include "hyperparameters_mvglm.stan" // same as hyperparameters.stan, but arrays of size M
}
transformed data {          
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);  // not used
  
  // declares poisson_max, hsM, idx_{global,local2,local4,mix,ool,noise}, 
  //   len_{global,local2,local4,mix,ool,noise}, {sqrt_,log_,sum_log_}y, 
  //   len_z_T, len_var_group, delta, is_continuous, pos, beta_smooth
  #include "tdata_mvglm.stan" 

  // defines hsM, idx_{global,local2,local4,mix,ool,noise}, 
  //   len_{global,local2,local4,mix,ool}, {sqrt_,log_,sum_log_}y, 
  //   len_z_T, len_var_group, delta, is_continuous, pos  
  #include "tdata2_mvglm.stan" 
}
parameters {
  // declares gamma_{nob,lob,upb}, z_beta, global, local{2,4}, mix, 
  //   ool, noise, aux_unscaled, z_b, z_T, rho, zeta, tau
  #include "parameters_mvglm.stan"
}
transformed parameters { 
  #include "tparameters_mvglm.stan" // defines aux, beta, b{_not_by_model}, theta_L
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 1.0, tau, scale, zeta, rho, z_T);
    b_not_by_model = make_b(z_b, theta_L, p, l);
    if (M == 1) b = b_not_by_model;
	else b = reorder_b(b_not_by_model, p, pmat, q1, q2, qmat, l, M);
  }
}
model {
  int aux_mark = 1; 
  #include "make_eta.stan" // defines eta
  if (t > 0) {
    #include "eta_add_Zb.stan"
  } 
  // Log-likelihoods
  for (m in 1:M) {
    vector[NM[m]] eta_tmp;	            // eta for just one submodel 
    eta_tmp = eta[idx[m,1]:idx[m,2]];   // eta for just one submodel 
    #include "eta_intercept_mvmer.stan"	// adds intercept or shifts eta
    if (family[m] == 8) {  // poisson-gamma mixture
	  #include "eta_add_noise_mvmer.stan"
    }    
    #include "mvmer_lp.stan" // increments target with long log-liks
    if (has_aux[m] == 1) aux_mark = aux_mark + 1;
  }
  // Log-priors	
  for (m in 1:M) { 
    if (has_aux[m] == 1) {
      aux_mark = sum(has_aux[1:m]);
      aux_lp(aux_unscaled[aux_mark], prior_dist_for_aux[m], prior_scale_for_aux[m], prior_df_for_aux[m])
    } 
  }
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  #include "gen_quantities_mvmer.stan"  // defines alpha, mean_PPD
}
