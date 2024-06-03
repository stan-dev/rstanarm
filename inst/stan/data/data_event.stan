  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus,
  //   5 = laplace, 6 = lasso
  int<lower=0,upper=6> e_prior_dist;
  int<lower=0,upper=2> e_prior_dist_for_intercept;

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> e_prior_dist_for_aux; // prior for basehaz params

  // data for event submodel
  real norm_const;            // constant shift for log baseline hazard
  int<lower=0> e_K;           // num. of predictors in event submodel
  int<lower=0> Npat;          // num. individuals (equal to l[id_var] - 1)
  int<lower=0> Nevents;       // num. events (ie. not censored)
  int<lower=0> qnodes;        // num. of nodes for Gauss-Kronrod quadrature
  int<lower=0> Npat_times_qnodes;
  int<lower=1,upper=3> basehaz_type;    // 1 = weibull, 2 = B-splines, 3 = piecewise
  int<lower=0> basehaz_df;              // df for baseline hazard
  int<lower=0,upper=1> e_has_intercept; // 1 = yes
  int<lower=0> nrow_e_Xq;     // num. rows in event predictor matrix at quad points
  matrix[e_K > 0 ? nrow_e_Xq : 0, e_K] e_Xq; // predictor matrix (event submodel) at qpts, centred
  vector[nrow_e_Xq] e_times;  // event times and unstandardised quadrature points
  matrix[nrow_e_Xq,basehaz_df] basehaz_X; // design matrix (basis terms) for baseline hazard
  vector[e_K] e_xbar;         // predictor means (event submodel)
  vector[Npat] e_weights;     // weights, set to zero if not used
  vector[Npat_times_qnodes] e_weights_rep; // repeated weights, set to zero if not used
  vector[Npat_times_qnodes] qwts; // GK quadrature weights with (b-a)/2 scaling
