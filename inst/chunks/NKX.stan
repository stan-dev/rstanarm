  // dimensions
  int<lower=0> N;  // number of observations
  int<lower=0> K;  // number of predictors
  
  // data
  vector[K] xbar;               // predictor means
  int<lower=0,upper=1> dense_X; // flag for dense vs. sparse
  matrix[N,K] X[dense_X];       // centered predictor matrix in the dense case

  // stuff for the sparse case
  int<lower=0> nnz_X;                    // number of non-zero elements in the implicit X matrix
  vector[nnz_X] w_X;                     // non-zero elements in the implicit X matrix
  int<lower=0> v_X[nnz_X];               // column indices for w_X
  int<lower=0> u_X[dense_X ? 0 : N + 1]; // where the non-zeros start in each row of X
