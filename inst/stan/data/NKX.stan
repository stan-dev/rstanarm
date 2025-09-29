  // dimensions
  int<lower=0> N;  // number of observations
  int<lower=0> K;  // number of predictors
  
  // data
  vector[K] xbar;               // predictor means
  int<lower=0,upper=1> dense_X; // flag for dense vs. sparse
  array[dense_X] matrix[N,K] X; // centered predictor matrix in the dense case

  // stuff for the sparse case
  int<lower=0> nnz_X;                        // number of non-zero elements in the implicit X matrix
  vector[nnz_X] w_X;                         // non-zero elements in the implicit X matrix
  array[nnz_X] int<lower=1, upper = K> v_X;  // column indices for w_X
  // where the non-zeros start in each row of X
  array[dense_X ? 0 : N + 1] int<lower=1, upper = rows(w_X) + 1> u_X; 
  
  // smooths
  int<lower=0> K_smooth;
  matrix[N,K_smooth] S;
  array[K_smooth] int<lower=1> smooth_map;
