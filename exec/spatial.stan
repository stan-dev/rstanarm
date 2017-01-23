// Spatial Autocorrelation Model (SEM)
// try demeaning predictors
// try non centered vs centered parameterization
data {
  int<lower=0> N;           // number of obs
  int<lower=0> K;           // number of predictors
  matrix[N,K] X;            // design matrix
  matrix<lower=0>[N,N] W;   // spatial weight matrix
  vector[N] y;              // response vector
  real<lower=1,upper=2> mod; // model type (SAR or SEM)
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
  
  if(mod == 1) {
    // target += multi_normal_prec_lpdf(y | inverse(weight_stuff) * X * beta, Sigma);
    target += -0.5 * quad_form(Sigma_inv, y - inverse(weight_stuff) * X * beta);
    target += -0.5 * log_determinant(Sigma_inv);
    target += cauchy_lpdf(rho | 0, 2);            // prior on spatial autocorrelation
    target += cauchy_lpdf(beta | 0, 3);           // priors on predictor parameters
    // target += cauchy_lpdf(sigma | 0, 3);          // prior on variation 
  }
}
generated quantities {
  real mean_PPD = 0;
  {
    vector[N] yrep;
    vector[N] eta;
    matrix[N,N] Sigma;
    eta = inverse(I - rho * W) * X * beta;
    Sigma = inverse(tcrossprod(I - rho * W));
    if(mod == 1)
      yrep = multi_normal_rng(eta, inverse(Sigma));
    mean_PPD = mean(yrep);
  }
}

