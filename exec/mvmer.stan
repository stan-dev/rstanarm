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
  // declares
  #include "dimensions_mvmer.stan"
  
	// declares
	#include "data_mvmer.stan" 

	// declares y_{,global_}prior_{mean,scale,df}{,_for_intercept,_for_aux}
  #include "hyperparameters_mvmer.stan"
}
transformed data {
  // declares
  #include "tdata_mvglm.stan" 
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
