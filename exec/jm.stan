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
  // declares
  #include "dimensions_mvmer.stan"
  
  // declares
  #include "data_mvmer.stan" 

  // declares e_prior_dist{_for_intercept,_for_aux}, Npat{_times_}qnodes, qwts, 
  //   basehaz_{type,df,X}, nrow_e_Xq, e_{K,Xq,times,d,xbar,weights,weights_rep}  
  #include "data_event.stan"

  // declares a_K, a_prior_dist, assoc, assoc_uses, has_assoc, {sum_}size_which_b, 
  //   which_b_zindex, {sum_}size_which_coef, which_coef_{zindex,xindex}, 
  //   {sum_}a_K_data, {sum_,sum_size_}which_interactions, y_Xq_{eta,eps,lag,auc,data},
  //   {nnz,w,v,u}_Zq_{eta,eps,lag,auc}, nrow_y_Xq, nrow_y_Xq_auc, 
  //   auc_qnodes, auc_qwts   
  #include "data_assoc.stan"
  
  // declares {e_,a_}{prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, 
  //   prior_{mean, scale, df}_for_aux, global_prior_{df,scale}}
  #include "hyperparameters_mvmer.stan"
  #include "hyperparameters_event.stan" 
  #include "hyperparameters_assoc.stan" 
}
transformed data {
  int<lower=0> e_hs = get_nvars_for_hs(e_prior_dist);
  int<lower=0> a_hs = get_nvars_for_hs(a_prior_dist);                 

  // declares poisson_max, hsM, idx_{global,local2,local4,mix,ool,noise}, 
  //   len_{global,local2,local4,mix,ool,noise}, {sqrt_,log_,sum_log_}y, 
  //   len_z_T, len_var_group, delta, is_continuous, pos, beta_smooth
  #include "tdata_mvmer.stan" 
}
parameters {
  // declares gamma_{nob,lob,upb}, z_beta, global, local{2,4}, mix, 
  //   ool, noise, aux_unscaled, z_b, z_T, rho, zeta, tau
  #include "parameters_mvmer.stan"
  
  // declares e_{gamma,z_beta,aux_unscaled,global,local,mix,ool}  
  #include "parameters_event.stan"
  
  // declares a_{z_beta,global,local,mix,ool}
  #include "parameters_assoc.stan"  
}
transformed parameters { 
  // declare parameters for event submodel
  vector[e_K] e_beta;               // log hazard ratios
  vector[a_K] a_beta;               // assoc params
  vector[basehaz_df] e_aux;         // basehaz params  
  
  // parameters for longitudinal submodels
  #include "tparameters_mvmer.stan" // defines aux, beta, b, theta_L
  
  // define parameters for event submodel
  e_beta = make_beta(e_z_beta, e_prior_dist, e_prior_mean, 
                     e_prior_scale, e_prior_df, e_global_prior_scale, 
                     e_global, e_local, e_ool, e_mix, rep_array(1.0, 0), 0);  
  a_beta = make_beta(a_z_beta, a_prior_dist, a_prior_mean, 
                     a_prior_scale, a_prior_df, a_global_prior_scale, 
                     a_global, a_local, a_ool, a_mix, rep_array(1.0, 0), 0);         
  e_aux  = make_basehaz_coef(e_aux_unscaled, e_prior_dist_for_aux,
                             e_prior_mean_for_aux, e_prior_scale_for_aux);
}
model {

  //---- Log likelihoods for longitudinal submodels
  
	#include "mvmer_lp.stan" // increments target with long log-liks
      
  //---- Log likelihood for event submodel (GK quadrature)
  {
	vector[nrow_e_Xq] e_eta_q; // eta for event submodel (at event and quad times)  
	
	// Event submodel: linear predictor at event and quad times
	if (e_K > 0) e_eta_q = e_Xq * e_beta;
	else e_eta_q = rep_vector(0.0, nrow_e_Xq);
	if (assoc == 1) { 
		// declares y_eta_q{_eps,_lag,_auc}, y_eta_qwide{_eps,_lag,_auc}, 
		//   y_q_wide{_eps,_lag,_auc}, mark{2,3}
		#include "assoc_evaluate.stan"
	}
			
	{
	// declares log_basehaz, log_{haz_q,haz_etimes,surv_etimes,event}
	#include "event_lp.stan" // increments target with event log-lik
	}
  }
    
  //---- Log priors

	#include "priors_mvmer.stan" // increments target with mvmer priors
	beta_lp(e_z_beta, e_prior_dist, e_prior_scale, e_prior_df, 
					e_global_prior_df, e_local, e_global, e_mix, e_ool);
	beta_lp(a_z_beta, a_prior_dist, a_prior_scale, a_prior_df, 
					a_global_prior_df, a_local, a_global, a_mix, a_ool);
	basehaz_lp(e_aux_unscaled, e_prior_dist_for_aux, 
						 e_prior_scale_for_aux, e_prior_df_for_aux);
	if (e_has_intercept == 1) 
		gamma_lp(e_gamma[1], e_prior_dist_for_intercept, e_prior_mean_for_intercept, 
						 e_prior_scale_for_intercept, e_prior_df_for_intercept);     
}
generated quantities {
  real e_alpha; // transformed intercept for event submodel 
    
  // declares and defines alpha, mean_PPD, cov matrix for lkj prior
  #include "gen_quantities_mvmer.stan"
    
  // norm_const is a constant shift in log baseline hazard
  if (e_has_intercept == 1) 
    e_alpha = e_gamma[1] + norm_const - 
      dot_product(e_xbar, e_beta) - dot_product(a_xbar, a_beta);
  else
    e_alpha = norm_const - 
      dot_product(e_xbar, e_beta) - dot_product(a_xbar, a_beta);
}
