  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> e_prior_dist;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  int<lower=0,upper=2> e_prior_dist_for_intercept;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> e_prior_dist_for_aux; // prior for basehaz params

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  int<lower=1,upper=3> basehaz_type;

  // dimensions
  int<lower=0> e_K;           // num. predictors in event submodel
  int<lower=0> basehaz_nvars; // num. aux parameters for baseline hazard
  int<lower=0> qnodes;        // num. nodes for GK quadrature
  int<lower=0> nevent;        // num. rows w/ an event (ie. not censored)
  int<lower=0> qevent;        // num. quadrature points for rows w/ an event
  int<lower=0> qlcens;        // num. quadrature points for rows w/ left censoring
  int<lower=0> qrcens;        // num. quadrature points for rows w/ right censoring
  int<lower=0> qicens;        // num. quadrature points for rows w/ interval censoring
  int<lower=0> qdelay;        // num. quadrature points for rows w/ delayed entry
  int<lower=0> idx_cpts[7,2]; // index for breaking cpts into epts,qpts_event,etc
  int<lower=0> len_cpts;

  // response and time variables
  vector[len_cpts] cpts;      // time at events and all quadrature points

  // predictor matrices
  matrix[len_cpts, e_K] e_x;  // predictor matrix
  vector[e_K] e_xbar;         // predictor means

  // spline basis for baseline hazard
  matrix[len_cpts, basehaz_nvars] basis_cpts;

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[qevent] qwts_event;
  vector[qlcens] qwts_lcens;
  vector[qrcens] qwts_rcens;
  vector[qicens] qwts_icenl;
  vector[qicens] qwts_icenu;
  vector[qdelay] qwts_delay;

  // constant shift for log baseline hazard
  real norm_const;

  // flags
  int<lower=0,upper=1> e_has_intercept; // basehaz requires intercept
