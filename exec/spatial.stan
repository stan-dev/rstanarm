// Spatial Autocorrelation Model (SEM)
// try demeaning predictors
// try non centered vs centered parameterization
data {
  int<lower=0> N;            // number of obs
  int<lower=0> K;            // number of predictors (excluding intercept)
  matrix[N,K] X;             // predictor matrix
  matrix<lower=0>[N,N] W;    // spatial weight matrix
  vector[N] y;               // response vector
  real<lower=1,upper=2> mod; // model type (1=SAR, 2=SEM)
  real<lower=0> shape1;               // alpha shape par for rho
  real<lower=0> shape2;               // beta shape par for rho
  int<lower=0,upper=1> has_intercept; // intercept included or not
  vector[K] xbar;                     // column means of X
}
transformed data {
  matrix[N,N] I;
  I = diag_matrix(rep_vector(1.0, N));
}
parameters {
  real<lower=0,upper=1> rho;     // note the support is (0,1) not (-1,1)
  vector[K] beta;                // parameter vector
  real gamma[has_intercept];     // intercept (demeaned)
}
model {
  matrix[N,N] weight_stuff;
  matrix[N,N] Sigma_inv;
  vector[N] eta;
  
  weight_stuff = I - rho * W;
  Sigma_inv = tcrossprod(weight_stuff);
  
  if(has_intercept == 1)
    eta = X * beta + gamma[1];
  else
    eta = X * beta;
    
  // model
  if(mod == 1) { // SAR
    // target += multi_normal_prec_lpdf(y | inverse(weight_stuff) * X * beta + , Sigma);
    target += -0.5 * quad_form(Sigma_inv, y - inverse(weight_stuff) * eta);
  }
  if(mod == 2) { // SEM
    target += multi_normal_prec_lpdf(y | eta, Sigma_inv);
  }
  // priors
  target += -0.5 * log_determinant(Sigma_inv);
  target += beta_lpdf(rho | shape1, shape2);       // prior on spatial autocorrelation
  // target += cauchy_lpdf(beta | 0, 3);           // priors on predictor parameters
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

