#include /pre/Columbia_copyright.stan
#include /pre/Brilleman_copyright.stan
#include /pre/license.stan

// Shared parameter joint model
functions {
#include /functions/common_functions.stan
#include /functions/bernoulli_likelihoods.stan
#include /functions/binomial_likelihoods.stan
#include /functions/continuous_likelihoods.stan
#include /functions/count_likelihoods.stan
#include /functions/mvmer_functions.stan
#include /functions/jm_functions.stan
}
data {
  // declares: M, has_aux, has_weights, resp_type, intercept_type,
  //   yNobs, yNeta, yK, t, p, l, q, len_theta_L, bN1, bK1, bK1_len
  //   bK1_idx, bN2, bK2, bK2_len, bK2_idx
#include /data/dimensions_mvmer.stan

  // declares: yInt{1,2,3}, yReal{1,2,3}, yX{1,2,3}, yXbar{1,2,3},
  //   family, link, y{1,2,3}_Z{1,2}, y{1,2,3}_Z{1,2}_id,
  //   y_prior_dist{_for_intercept,_for_aux,_for_cov}, prior_PD
#include /data/data_mvmer.stan

  // declares: e_prior_dist{_for_intercept,_for_aux},
  //   Npat, Nevents, qnodes, Npat_times_qnodes, qwts,
  //   basehaz_{type,df,X}, nrow_e_Xq, e_has_intercept, nrow_e_Xq,
  //   e_{K,Xq,times,xbar,weights,weights_rep}
#include /data/data_event.stan

  // declares: a_{K,xbar}, a_prior_dist, assoc, assoc_uses, has_assoc,
  //   {sum_}size_which_b, which_b_zindex, {sum_}size_which_coef,
  //   which_coef_{zindex,xindex}, a_K_data, y_Xq_{eta,eps,lag,auc,data},
  //   {sum_,sum_size_}which_interactions, idx_q,
  //   nrow_y_Xq{_auc}, auc_{qnodes,qwts}, has_grp, grp_assoc, grp_idx,
  //   y{1,2,3}_xq_{eta,eps,auc}, y{1,2,3}_z{1,2}q_{eta,eps,auc},
  //   y{1,2,3}_z{1,2}q_id_{eta,eps,auc}
#include /data/data_assoc.stan

  // declares: e_prior_{mean,scale,df}{_for_intercept,for_aux},
  //   e_global_prior_{scale,df}
#include /data/hyperparameters_mvmer.stan
#include /data/hyperparameters_event.stan
  // declares: a_prior_{mean,scale,df}, a_global_prior_{scale,df}
#include /data/hyperparameters_assoc.stan
}
transformed data {
  int<lower=0> e_hs = get_nvars_for_hs(e_prior_dist);
  int<lower=0> a_hs = get_nvars_for_hs(a_prior_dist);

  // declares: yHs{1,2,3}, len_{z_T,var_group,rho}, pos, delta,
  //   bCov{1,2}_idx, {sqrt,log,sum_log}_y{1,2,3},
#include /tdata/tdata_mvmer.stan
}
parameters {
  // declares: yGamma{1,2,3}, z_yBeta{1,2,3}, z_b, z_T, rho,
  //   zeta, tau, bSd{1,2}, z_bMat{1,2}, bCholesky{1,2},
  //   yAux{1,2,3}_unscaled, yGlobal{1,2,3}, yLocal{1,2,3},
  //   yOol{1,2,3}, yMix{1,2,3}
#include /parameters/parameters_mvmer.stan

  // declares e_{gamma,z_beta,aux_unscaled,global,local,mix,ool}
#include /parameters/parameters_event.stan

  // declares a_{z_beta,global,local,mix,ool}
#include /parameters/parameters_assoc.stan
}
transformed parameters {
  vector[e_K] e_beta;               // log hazard ratios
  vector[a_K] a_beta;               // assoc params
  vector[basehaz_df] e_aux;         // basehaz params

  //---- Parameters for longitudinal submodels

  // declares and defines: yBeta{1,2,3}, yAux{1,2,3}, yAuxMaximum,
  //   theta_L, bMat{1,2}
#include /tparameters/tparameters_mvmer.stan
  //---- Parameters for event submodel
  e_beta = make_beta(e_z_beta, e_prior_dist, e_prior_mean,
                     e_prior_scale, e_prior_df, e_global_prior_scale,
                     e_global, e_local, e_ool, e_mix, rep_array(1.0, 0), 0,
                     e_slab_scale, e_caux);
  a_beta = make_beta(a_z_beta, a_prior_dist, a_prior_mean,
                     a_prior_scale, a_prior_df, a_global_prior_scale,
                     a_global, a_local, a_ool, a_mix, rep_array(1.0, 0), 0,
                     a_slab_scale, a_caux);
  e_aux  = make_basehaz_coef(e_aux_unscaled, e_prior_dist_for_aux,
                             e_prior_mean_for_aux, e_prior_scale_for_aux);
}
model {

  //---- Log likelihoods for longitudinal submodels
#include /model/mvmer_lp.stan

  //---- Log likelihood for event submodel (GK quadrature)
  {
    vector[nrow_e_Xq] e_eta_q; // eta for event submodel (at event and quad times)

    // Event submodel: linear predictor at event and quad times
    if (e_K > 0) e_eta_q = e_Xq * e_beta;
    else e_eta_q = rep_vector(0.0, nrow_e_Xq);
    if (assoc == 1) {
        // declares y_eta_q{_eps,_lag,_auc}, y_eta_qwide{_eps,_lag,_auc},
        //   y_q_wide{_eps,_lag,_auc}, mark{2,3}
#include /model/assoc_evaluate.stan
    }

    {
    // declares log_basehaz, log_{haz_q,haz_etimes,surv_etimes,event}
    // increments target with event log-lik
#include /model/event_lp.stan
    }
  }

  //---- Log priors
  // increments target with mvmer priors
#include /model/priors_mvmer.stan
    beta_lp(e_z_beta, e_prior_dist, e_prior_scale, e_prior_df,
                    e_global_prior_df, e_local, e_global, e_mix, e_ool,
                    e_slab_df, e_caux);
    beta_lp(a_z_beta, a_prior_dist, a_prior_scale, a_prior_df,
                    a_global_prior_df, a_local, a_global, a_mix, a_ool,
                    a_slab_df, a_caux);
    basehaz_lp(e_aux_unscaled, e_prior_dist_for_aux,
                         e_prior_scale_for_aux, e_prior_df_for_aux);
    if (e_has_intercept == 1)
        gamma_lp(e_gamma[1], e_prior_dist_for_intercept, e_prior_mean_for_intercept,
                         e_prior_scale_for_intercept, e_prior_df_for_intercept);
}
generated quantities {
  real e_alpha; // transformed intercept for event submodel

  // declares and defines: mean_PPD, yAlpha{1,2,3}, b{1,2}, bCov{1,2}
#include /gqs/gen_quantities_mvmer.stan

  // norm_const is a constant shift in log baseline hazard
  if (e_has_intercept == 1)
    e_alpha = e_gamma[1] + norm_const -
      dot_product(e_xbar, e_beta) - dot_product(a_xbar, a_beta);
  else
    e_alpha = norm_const -
      dot_product(e_xbar, e_beta) - dot_product(a_xbar, a_beta);
}
