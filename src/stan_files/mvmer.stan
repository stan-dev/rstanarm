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
#include /functions/jm_functions.stan
}
data {
  // declares
#include /data/dimensions_mvmer.stan

  // declares
#include /data/data_mvmer.stan
  
  // declares y_{,global_}prior_{mean,scale,df}{,_for_intercept,_for_aux}
#include /data/hyperparameters_mvmer.stan
}
transformed data {
  // declares
#include /tdata/tdata_mvmer.stan
}
parameters {
  // declares
#include /parameters/parameters_mvmer.stan
}
transformed parameters { 
  // declares
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
  // declares and defines alpha, mean_PPD, cov matrix for lkj prior
#include /gqs/gen_quantities_mvmer.stan
}
