  int<lower=0> num_non_zero;    // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;       // non-zero elements in the implicit Z matrix
  int<lower=0> v[num_non_zero]; // column indices for w
  int<lower=0> u[(N+1)*(t>0)];  // where the non-zeros start in each row
