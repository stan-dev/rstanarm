  // dimensions
  int<lower=0> N;  // number of observations
  int<lower=0> K;  // number of predictors
  
  // data
  vector[K] xbar;  // predictor means
  matrix[N,K] X;  // centered predictor matrix
