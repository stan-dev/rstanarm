  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes
  int<lower=0,upper=1> compute_mean_PPD; // 1 = yes
  
  // intercept
  int<lower=0,upper=1> has_intercept;  // 1 = yes
  
  // link function from location to linear predictor 
  int<lower=1> link;  // interpretation varies by .stan file
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> prior_dist;
  int<lower=0,upper=2> prior_dist_for_intercept;
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux;
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_smooth;
