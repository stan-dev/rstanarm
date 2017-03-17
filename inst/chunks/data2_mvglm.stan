  // combined outcome vectors
  vector[N_real] y_real;      // reals                
  int<lower=0> y_int[N_int];  // integers

  // num. binomial trials, set to 0 for obs. that are not binomial
  int<lower=0> trials[N];

  // auxiliary parameters
  int<lower=0,upper=1> has_aux[M];  // 1 = Yes
