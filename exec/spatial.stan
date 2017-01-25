// Spatial Autocorrelation Model (SEM)
// try demeaning predictors
// try non centered vs centered parameterization
data {
  int<lower=0> N;            // number of obs
  int<lower=0> K;            // number of predictors
  matrix[N,K] X;             // design matrix
  matrix<lower=0>[N,N] W;    // spatial weight matrix
  vector[N] y;               // response vector
  real<lower=1,upper=2> mod; // model type (SAR or SEM)
  real<lower=0> shape1;               // alpha shape par for rho
  real<lower=0> shape2;               // beta shape par for rho
}
transformed data {
  matrix[N,N] I;
  I = diag_matrix(rep_vector(1.0, N));
}
parameters {
  real<lower=0,upper=1> rho;     // note the support is (0,1) not (-1,1)
  vector[K] beta;                // parameter vector
  // real<lower=0> sigma;        // variance
}
model {
  matrix[N,N] weight_stuff;
  matrix[N,N] Sigma_inv;
  weight_stuff = I - rho * W;
  Sigma_inv = tcrossprod(weight_stuff);
  
  // model
  if(mod == 1) { // SAR
    // target += multi_normal_prec_lpdf(y | inverse(weight_stuff) * X * beta, Sigma);
    target += -0.5 * quad_form(Sigma_inv, y - inverse(weight_stuff) * X * beta);
  }
  if(mod == 2) { // SEM
    target += multi_normal_prec_lpdf(y | X * beta, Sigma_inv);
  }
  // priors
  target += -0.5 * log_determinant(Sigma_inv);
  target += beta_lpdf(rho | shape1, shape2);       // prior on spatial autocorrelation
  // target += cauchy_lpdf(beta | 0, 3);           // priors on predictor parameters
}
generated quantities {
  real mean_PPD = 0;
  {
    vector[N] yrep;
    vector[N] eta;
    matrix[N,N] Sigma;
    if(mod == 1)
      eta = inverse(I - rho * W) * X * beta;
    else if(mod == 2)
      eta = X * beta;
    Sigma = inverse(tcrossprod(I - rho * W));
    yrep = multi_normal_rng(eta, inverse(Sigma));
    mean_PPD = mean(yrep);
  }
}

