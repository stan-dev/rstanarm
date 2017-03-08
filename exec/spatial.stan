// Spatial Autocorrelation Model (SEM)
// try demeaning predictors
// try non centered vs centered parameterization
data {
  int<lower=0> N;                     // number of obs
  int<lower=0> K;                     // number of predictors (excluding intercept)
  matrix[N,K] X;                      // predictor matrix
  matrix<lower=0>[N,N] W;             // spatial weight matrix
  vector[N] y;                        // response vector
  real<lower=1,upper=2> mod;          // model type (1=SAR, 2=SEM)
  real<lower=0> shape1;               // alpha shape par for rho
  real<lower=0> shape2;               // beta shape par for rho
  int<lower=0,upper=1> has_intercept; // intercept included or not
  vector[K] xbar;                     // column means of X
  // prior info for intercept
  int<lower=0,upper=2> prior_dist_for_intercept;
  real prior_mean_for_intercept;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_df_for_intercept;
}
transformed data {
  matrix[N,N] I;
  I = diag_matrix(rep_vector(1.0, N));
}
parameters {
  real<lower=0,upper=1> rho;          // note the support is (0,1) not (-1,1)
  vector[K] beta;                     // parameter vector
  real gamma[has_intercept];          // intercept (demeaned)
  real<lower=0> sigma;
}
model {
  matrix[N,N] weight_stuff;
  vector[N] eta;
  vector[N] half;
  
  weight_stuff = I - rho * W;
  // intercept present
  if(has_intercept == 1)
    eta = X * beta + gamma[1];
  else
    eta = X * beta;
  // components of likelihood
  if(mod == 1)  // SAR
    half = weight_stuff * (y - mdivide_left(weight_stuff, eta));
  if(mod == 2)  // SEM
    half = weight_stuff * (y - eta);
  // likelihood
  target += log_determinant(weight_stuff) - N * log(sigma) - 0.5 * dot_self(half);
  // prior on spatial autocorrelation
  target += beta_lpdf(rho | shape1, shape2);
  // prior on intercept
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1)  // normal
      target += normal_lpdf(gamma[1] | prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2)  // student_t
      target += student_t_lpdf(gamma[12] | prior_df_for_intercept, prior_mean_for_intercept, 
                               prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD = 0;
  // scale intercept if included
    
  // calculate mean of posterior predictive distribution
  {
    vector[N] yrep;
    vector[N] eta;
    vector[N] eta2;
    matrix[N,N] Sigma;
    // fix intercept if included
    if (has_intercept == 1) {
      alpha[1] = gamma[1] - dot_product(beta, xbar);
      eta = X * beta + alpha[1];
    }
    else {
      eta = X * beta;
    }
    // define model specific linear predictor
    if(mod == 1)
      eta2 = inverse(I - rho * W) * eta;
    else if(mod == 2)
      eta2 = eta;
    Sigma = inverse(tcrossprod(I - rho * W));
    yrep = multi_normal_rng(eta2, inverse(Sigma));
    mean_PPD = mean(yrep);
  }
}

