  // combined outcome vectors
  vector[N_real] y_real;      // reals                
  int<lower=0> y_int[N_int];  // integers

  // auxiliary parameters
  int<lower=0,upper=1> has_aux[M];  // 1 = Yes
