  int<lower=0> num_non_zero;  // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;     // non-zero elements in the implicit Z matrix
  int<lower=0, upper=q-1> v[num_non_zero];               // column indices for w
  int<lower=0, upper=rows(w) + 1> u[t > 0 ? N + 1 : 0];  // where the non-zeros start in each row
  int<lower=0,upper=1> special_case;                     // is the only term (1|group)
