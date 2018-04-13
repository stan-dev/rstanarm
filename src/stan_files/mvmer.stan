#include /pre/Columbia_copyright.stan
#include /pre/Brilleman_copyright.stan
#include /pre/license.stan

// Multivariate GLM with correlated group-specific terms
functions {
#include /functions/common_functions.stan
#include /functions/bernoulli_likelihoods.stan
#include /functions/binomial_likelihoods.stan
#include /functions/continuous_likelihoods.stan
#include /functions/count_likelihoods.stan
#include /functions/mvmer_functions.stan
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
  
  // declares: y_prior_{mean,scale,df}{1,2,3,_for_intercept,_for_aux}, 
  //   y_global_prior_{df,scale}, len_{concentration,regularization},
  //   b_prior_{shape,scale,concentration,regularization},
  //   b{1,2}_prior_{scale,df,regularization}
#include /data/hyperparameters_mvmer.stan
}
transformed data {
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
}
transformed parameters { 
  // declares and defines: yBeta{1,2,3}, yAux{1,2,3}, yAuxMaximum, 
  //   theta_L, bMat{1,2}
#include /tparameters/tparameters_mvmer.stan
}
model {
  // Log likelihoods
  // increments target with mvmer log liks
#include /model/mvmer_lp.stan

  // Log priors
  // increments target with mvmer priors
#include /model/priors_mvmer.stan
}
generated quantities {
  // declares and defines: mean_PPD, yAlpha{1,2,3}, b{1,2}, bCov{1,2}
#include /gqs/gen_quantities_mvmer.stan
}
