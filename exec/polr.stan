# GLM for an ordinal outcome with coherent priors
data {
  int<lower=2> J;                # number of outcome categories, which typically is > 2
  int<lower=1> N;                # number of observations
  int<lower=0> K;                # number of predictors (excluding a constant)
  matrix[N,K]  X;                # centered predictor matrix
  vector[K] xbar;                # means of the predictors
  vector<lower=0>[K] s_X;        # standard deviations of the predictors
  int<lower=1,upper=J> y[N];     # outcome
  int<lower=0,upper=1> prior_PD; # flag indicating whether to draw from the prior predictive
  
  # link function from location to linear predictor
  int<lower=1,upper=5> link;

  # weights
  int<lower=0,upper=1> has_weights; # 0 = No, 1 = Yes
  vector[N * has_weights] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No, 1 = Yes
  vector[N * has_offset] offset;
  
  # prior family (zero indicates no prior!!!)
  int<lower=0,upper=1> prior_dist;  # 1 = Ben
  
  # hyperparameter values
  real<lower=0> shape;
  vector<lower=0>[J] prior_counts;
  
  int<lower=0,upper=1> do_residuals;
}
transformed data {
  real<lower=shape> shapephalf;
  real<lower=0> half_K;
  int<lower=0,upper=1> is_constant;
  matrix[K,K] middle;
  shapephalf <- shape + 0.5;
  half_K <- 0.5 * K;
  is_constant <- 1;
  for (j in 1:J) if (prior_counts[j] != 1) is_constant <- 0;
  middle <- xbar * transpose(xbar);
}
parameters {
  simplex[J] pi;
  row_vector[K] z_beta;
  cholesky_factor_corr[K] L[prior_dist == 1];
  real<lower=0,upper=1>  R2[prior_dist == 1];
}
transformed parameters {
  real Delta_y;
  vector[K] beta;
  vector[J-1] cutpoints;
  if (prior_dist == 1) {
    Delta_y <- inv(sqrt(1 - R2[1]));
    if (K > 1)
      beta <- transpose(mdivide_right_tri_low(z_beta, L[1])) *
              sqrt(R2[1] / dot_self(z_beta)) ./ s_X * Delta_y;
    else beta[1] <- sqrt(R2[1]) / s_X[1] * Delta_y;
  }
  else { // prior_dist == 0
    beta <- transpose(z_beta);
    Delta_y <- sqrt(quad_form(middle, beta) + 1);
  }
  cutpoints <- make_cutpoints(pi, Delta_y, link);
}
model {
  vector[N] eta;
  if (K > 0) eta <- X * beta;
  else eta <- rep_vector(0.0, N);
  if (has_offset == 1) eta <- eta + offset;
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    increment_log_prob(pw_polr(y, eta, cutpoints, link));
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    increment_log_prob(dot_product(weights, pw_polr(y, eta, cutpoints, link)));
  }
  
  if (prior_dist == 1) {
    z_beta ~ normal(0, 1);
    if (K > 1) L[1] ~ lkj_corr_cholesky(shapephalf);
    R2[1] ~ beta(half_K, shape);
    if (is_constant == 0) pi ~ dirichlet(prior_counts);
  }
  /* else prior_dist is 0 and nothing is added */
}
generated quantities {
  vector[J-1] zeta;
  vector[(J > 2) * (J - 1) + 1] mean_PPD;
  vector[N * do_residuals] residuals;
  zeta <- cutpoints + dot_product(xbar, beta);
  if (J == 2) zeta <- -zeta;
  mean_PPD <- rep_vector(0,rows(mean_PPD));
  {
    vector[N] eta;
    if (K > 0) eta <- X * beta;
    else eta <- rep_vector(0.0, N);
    if (has_offset == 1) eta <- eta + offset;
    for (n in 1:N) {
      vector[J] theta;
      int y_tilde;
      real previous;
      real ystar;
      theta[1] <- CDF_polr(cutpoints[1] - eta[n], link);
      previous <- theta[1];
      for (j in 2:(J-1)) {
        real current;
        current <- CDF_polr(cutpoints[j] - eta[n], link);
        theta[j] <- current - previous;
        previous <- current;
      }
      theta[J] <- 1 - previous;
      if (J == 2) {
        mean_PPD[1] <- mean_PPD[1] + bernoulli_rng(theta[J]);
      }
      else {
        y_tilde <- categorical_rng(theta);
        mean_PPD[y_tilde] <- mean_PPD[y_tilde] + 1;
      }
      
      if (do_residuals) {
        if (y[n] == 1)
          ystar <- draw_ystar_rng(negative_infinity(), cutpoints[1], eta[n], link);
        else if (y[n] == J)
          ystar <- draw_ystar_rng(cutpoints[J - 1], positive_infinity(), eta[n], link);
        else ystar <- draw_ystar_rng(cutpoints[y[n] - 1], cutpoints[y[n]], eta[n], link);
        residuals[n] <- ystar - eta[n];
      }
    }
    mean_PPD <- mean_PPD / N;
  }
}
