  int<lower=1> M;       // num. of long. submodels
  int<lower=0> NM[M];   // num. of obs. in each submodel
  int<lower=0> N_real;  // num. of obs. across all submodels with real outcomes
  int<lower=0> N_int;   // num. of obs. across all submodels with integer outcomes
  int<lower=0> N01[M,2];                   // num. of bernoulli 0/1 observations in each submodel
  int<lower=0,upper=N> idx[M,2];           // indices of first and last obs. for each submodel
  int<lower=0,upper=N_real> idx_real[M,2]; // indices of first and last real obs. for each submodel
  int<lower=0,upper=N_int>  idx_int [M,2]; // indices of first and last integer obs. for each submodel
  int<lower=0> KM[M];                      // num. of predictors in each submodel
  int<lower=0,upper=K> idx_K[M,2];         // indices of first and last betas for each submodel
  int<lower=0,upper=M> sum_has_intercept;     // num. submodels w/ intercept
  int<lower=0,upper=M> sum_has_intercept_nob; // num. submodels w/ unbounded intercept
  int<lower=0,upper=M> sum_has_intercept_lob; // num. submodels w/ lower bounded intercept
  int<lower=0,upper=M> sum_has_intercept_upb; // num. submodels w/ upper bounded intercept
  int<lower=0,upper=M> sum_has_aux;           // num. submodels w/ auxiliary term
