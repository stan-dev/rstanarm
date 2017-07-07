  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes
  
  // intercepts
  int<lower=0,upper=1> has_intercept[M];      // 1 = yes
  int<lower=0,upper=1> has_intercept_nob[M];  // unbounded
  int<lower=0,upper=1> has_intercept_lob[M];  // lower bound at 0
  int<lower=0,upper=1> has_intercept_upb[M];  // upper bound at 0

  // family
  int<lower=1> family[M];  
 
  // link function, varies by family 
  int<lower=1> link[M];    

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> prior_dist[M];
  int<lower=0,upper=2> prior_dist_for_intercept[M];  
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux[M];
