#include "Columbia_copyright.stan"
#include "Brilleman_copyright.stan"
#include "license.stan" // GPL3+

// Multivariate GLM with correlated group-specific terms
functions {
  #include "common_functions.stan"
  #include "bernoulli_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "continuous_likelihoods.stan"
  #include "count_likelihoods.stan"
  #include "jm_functions.stan"
}
data {
  // declares: M, has_aux, has_weights, resp_type, intercept_type,
	//   yNobs, yNeta, yK, t, p, l, q, len_theta_L, bN1, bK1, bK1_len
	//   bK1_idx, bN2, bK2, bK2_len, bK2_idx
  #include "dimensions_mvmer.stan"

  // declares: yInt{1,2,3}, yReal{1,2,3}, yX{1,2,3}, yXbar{1,2,3},
	//   family, link, y{1,2,3}_Z{1,2}, y{1,2,3}_Z{1,2}_id,
	//   y_prior_dist{_for_intercept,_for_aux,_for_cov}, prior_PD
  #include "data_mvmer.stan"
  
	// declares: y_prior_{mean,scale,df}{1,2,3,_for_intercept,_for_aux}, 
	//   y_global_prior_{df,scale}, len_{concentration,regularization},
	//   b_prior_{shape,scale,concentration,regularization},
	//   b{1,2}_prior_{scale,df,regularization}
  #include "hyperparameters_mvmer.stan"
}
transformed data {
  // declares
  #include "tdata_mvmer.stan" 
}
parameters {
  // declares
  #include "parameters_mvmer.stan"
}
transformed parameters { 
  // declares
  #include "tparameters_mvmer.stan"
}
model {
  // Log likelihoods
  #include "mvmer_lp.stan" // increments target with mvmer log liks

  // Log priors
  #include "priors_mvmer.stan" // increments target with mvmer priors
}
generated quantities {
  // declares and defines alpha, mean_PPD, cov matrix for lkj prior
  #include "gen_quantities_mvmer.stan"
}
