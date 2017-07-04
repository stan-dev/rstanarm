#include "Columbia_copyright.stan"
#include "Brilleman_copyright.stan"
#include "license.stan" // GPL3+

// Shared parameter joint model
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
  int<lower=1,upper=t> t_i;    // index of grouping factor corresponding to patient-level

  // declares e_prior_dist{_for_intercept,_for_aux}, Npat{_times_}quadnodes, quadweight, 
  //   basehaz_{type,df,X}, nrow_{y,e}_Xq, e_{K,Xq,times,d,xbar,weights,weights_rep}  
  #include "data_event.stan"

  // declares a_K, a_prior_dist, assoc, assoc_uses, has_assoc, {sum_}size_which_b, 
  //   which_b_zindex, {sum_}size_which_coef, which_coef_{zindex,xindex}, 
  //   {sum_}a_K_data, {sum_,sum_size_}which_interactions, y_Xq_{eta,eps,lag,auc,data},
  //   {nnz,w,v,u}_Zq_{eta,eps,lag,auc}, nrow_y_Xq_auc, auc_quadnodes, auc_quadweights   
  #include "data_assoc.stan"
  
  // declares {e_,a_}{prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, 
  //   prior_{mean, scale, df}_for_aux, global_prior_{df,scale}}
  #include "hyperparameters_mvglm.stan" // same as hyperparameters.stan, but arrays of size M
  #include "hyperparameters_event.stan" 
  #include "hyperparameters_assoc.stan" 
}
transformed data {
  int<lower=0> e_hs = get_nvars_for_hs(e_prior_dist);                 
  int<lower=0> a_hs = get_nvars_for_hs(a_prior_dist);                 
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
  // declares e_{gamma,z_beta,aux_unscaled,global,local,mix,ool}  
  #include "parameters_event.stan"
  // declares a_{z_beta,global,local,mix,ool}
  #include "parameters_assoc.stan"  
}
transformed parameters { 
  // parameters for event submodel
  vector[e_K] e_beta;               // log hazard ratios
  vector[a_K] a_beta;               // assoc params
  vector[basehaz_df] e_aux;         // basehaz params  
  #include "tparameters_mvglm.stan" // defines aux, beta, b{_not_by_model}, theta_L
  e_beta = generate_beta(e_z_beta, e_prior_dist, e_prior_mean, 
                         e_prior_scale, e_prior_df, e_global, e_local,
                         e_global_prior_scale, e_ool, e_mix);  
  a_beta = generate_beta(a_z_beta, a_prior_dist, a_prior_mean, 
                         a_prior_scale, a_prior_df, a_global, a_local,
                         a_global_prior_scale, a_ool, a_mix);         
  e_aux  = generate_aux(e_aux_unscaled, e_prior_dist_for_aux,
                        e_prior_mean_for_aux, e_prior_scale_for_aux);
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 1.0, tau, scale, zeta, rho, z_T);
    b_not_by_model = make_b(z_b, theta_L, p, l);
    if (M == 1) b = b_not_by_model;
	else b = reorder_b(b_not_by_model, p, pmat, q1, q2, qmat, l, M);
  }
}
model {
  vector[nrow_e_Xq] e_eta_q; // eta for event submodel (at event and quad times)  

  //---- Log-lik for longitudinal submodels
  
  int aux_mark = 1; 
  #include "make_eta.stan" // defines eta
  if (t > 0) {
    #include "eta_add_Zb.stan"
  } 
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
  
  //----- Log-lik for event submodel (GK quadrature)
  
  // Event submodel: linear predictor at event and quad times
  if (e_K > 0) e_eta_q = e_Xq * e_beta;
  else e_eta_q = rep_vector(0.0, nrow_e_Xq);
  if (e_has_intercept == 1) e_eta_q = e_eta_q + e_gamma[1];
  else e_eta_q = e_eta_q + dot_product(e_xbar, e_beta);     
  if (assoc == 1) { 
    // declares y_eta_q{_eps,_lag,_auc}, y_eta_qwide{_eps,_lag,_auc}, 
	  //   y_q_wide{_eps,_lag,_auc}, mark{2,3}
    #include "assoc_definitions.stan"  
    #include "assoc_prepwork.stan"
    #include "assoc_evaluate.stan"
  }
  { 
    // declares log_basehaz, ll_{haz_q,haz_eventtime,surv_eventtime,event}
	  #include "event_lp.stan" // increments target with event log-lik
  }
  
  //----- Log-priors

  // increment target with priors for aux params 
  for (m in 1:M) { 
    if (has_aux[m] == 1) {
      aux_mark = sum(has_aux[1:m]);
      aux_lp(aux_unscaled[aux_mark], prior_dist_for_aux[m], prior_scale_for_aux[m], prior_df_for_aux[m])
    } 
  }
  
  // increment target with priors for betas and gammas 
  #include "priors_mvglm.stan"  
  beta_lp(e_z_beta, e_prior_dist, e_prior_scale, e_prior_df, 
          e_global_prior_df, e_local, e_global, e_mix, e_ool)
  beta_lp(a_z_beta, a_prior_dist, a_prior_scale, a_prior_df, 
          a_global_prior_df, a_local, a_global, a_mix, a_ool)
  if (e_has_intercept == 1) 
    gamma_lp(e_gamma[1], e_prior_dist_for_intercept, e_prior_mean_for_intercept, 
             e_prior_scale_for_intercept, e_prior_df_for_intercept);  
   
  // increment target with priors for basehaz params
  basehaz_lp(e_aux_unscaled, e_prior_dist_for_aux, 
             e_prior_scale_for_aux, e_prior_df_for_aux)
  
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  #include "gen_quantities_mvmer.stan"  // defines alpha, mean_PPD
}
