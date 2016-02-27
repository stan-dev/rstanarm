functions {
  /** 
  * log-likelihood vector for a multinomial probit model
  *
  * The basic idea here is to introduce p - 1 unknown multinormal errors
  * that correspond to all the choices that correspond to all the choices
  * that unit i did not make. Then we can figure out a lower bound to the
  * error term that corresponds to the choice that unit i made, and integrate
  * the univariate normal distribution from that lower bound to infinity to
  * produce the likelihood of unit i making that choice.
  *
  * Everything is relative to the baseline choice to overcome location
  * shifts but the error correlation matrix is p x p. The error variances
  * are proportions of p to overcome rescalings. We then have to take
  * the difference in errors relative to the first.
  *
  * @param alpha Vector of intercepts
  * @param beta Matrix of regression coefficients on unit-specific predictors
  * @param gamma Vector of regression coefficients on alternative-specific predictors
  * @param L Cholesky factor of the full error correlation matrix
  * @param z_epsilon Array of vectors of primitives of the p - 1 errors
  * @param pm1 Integer equal to p - 1
  * @param X Array of vectors of centered unit-specific predictors
  * @param Z Two-dimensional array of vectors of centered alternative-specific 
  *          predictors
  * @param y Integer array of observed choices
  * @param include Two-dimensional integer array of indices
  * @return A scalar log-likelihood
  */
  
  real ll_mnp(vector alpha, matrix beta, vector gamma,
              matrix L, vector z_epsilon,
              int pm1, vector[] X, vector[,] Z, int[] y, int[,] include) {
                       
    matrix[pm1,pm1] Lp[pm1 + 1]; // array of Cholesky factors
    vector[size(y)] lowers;      // lower bounds to the normal CCDF
    vector[size(y)] sigma;       // standard deviations for the normal CCDF
    int p;                       // number of choices
    real sqrt_p;
    vector[1] zero;
    vector[pm1 + 1] diff_sds;    // standard deviation in utility relative to choice 1
    int pos;
    p <- pm1 + 1;
    sqrt_p <- sqrt(p);
    zero[1] <- 0;
    {
      matrix[p,p] Lambda;         // full covariance matrix among the errors
      Lambda <- multiply_lower_tri_self_transpose(L);
      for (i in 1:pm1) Lp[i] <- cholesky_decompose(Lambda[include[i],include[i]]);
      Lp[p] <- L[include[p],include[p]];
      for (i in 1:pm1) diff_sds[i] <- sqrt(2 - 2 * Lambda[i+1,1]);
    }

    pos <- 1;
    for (i in 1:size(y)) {      // specify the lower bound for each chosen error
      vector[pm1] v;
      vector[pm1] epsilon;
      int y_i;
      y_i <- y[i];
      if (rows(alpha) > 0) v <- alpha;
      else v <- rep_vector(0, pm1);
      if (rows(beta) > 0) v <- v + beta * X[i];
      for (j in 1:rows(gamma)) v <- v + gamma[j] * Z[i,j];
      epsilon <- Lp[y_i] * segment(z_epsilon, pos, pm1);
      pos <- pos + pm1;
      if (y_i == 1) lowers[i] <- max(v + epsilon);
      else lowers[i] <- max(append_row(zero, v[include[y_i]]) + epsilon) + v[y_i];
      sigma[i] <- diff_sds[y_i];
    }
    return normal_ccdf_log(lowers, 0.0, sigma);
  }
  
}
data {
  int<lower=1> N;                      // number of observations
  int<lower=0> K;                      // number of unit-specific predictors
  vector[K] X[N];                      // centered  unit-specific predictors
  int<lower=3> p;                      // number of choices
  int<lower=1,upper=p> y[N];           // observed choice outcomes
  int<lower=0> q;                      // number of alternative-specific predictors
  vector[p - 1] Z[N, q];               // centered alternative-specific predictors

  vector<lower=0>[p] prior_counts;     // hyperparameters
  real<lower=0> eta;
  int<lower=3> sigma_vech_size;
}
transformed data {
  matrix[p,p-1] Tbc_t;                 // transformation matrix relative to baseline
  int<lower=1,upper=p> include[p,p-1];
  int<lower=0,upper=1> flat;
  int<lower=2> pm1;
  real<lower=0> sqrt_p;
  sqrt_p <- sqrt(p);
  pm1 <- p - 1;
  flat <- 1;
  for (i in 1:p) {
    if (prior_counts[i] != 1) flat <- 0;
    if (i > 1) for (j in 1:(i-1)) include[i,j] <- j;
    if (i < p) for (j in (i+1):p) include[i,j-1] <- j;
  }
  Tbc_t <- append_row(rep_row_vector(-1, pm1), 
                      diag_matrix(rep_vector(1, pm1)));
}
parameters {
  vector[pm1] alpha;                   // intercepts
  matrix[pm1,K] beta;                  // coefficients on unit-specific predictors
  vector[q] gamma;                     // coefficients on alternative-specific predictors
  
  cholesky_factor_corr[p] L;           // Cholesky factor of full error correlation matrix
  vector[pm1 * N] z_epsilon;           // primitives of the errors
}
model {
  increment_log_prob(ll_mnp(alpha, beta, gamma, 
                            L, z_epsilon,
                            pm1, X, Z, y, include));
  // priors
  // prior on alpha
  // prior on beta
  L ~ lkj_corr_cholesky(eta);
  z_epsilon ~ normal(0,1);
}
generated quantities {
  matrix[pm1,pm1] Lambda;
  vector[p] mean_PPD;                  // average posterior predictive distribution
  mean_PPD <- rep_vector(0, p);
  {
    matrix[pm1,pm1] relative_L;
    vector[pm1] zeros;
    Lambda <- quad_form(multiply_lower_tri_self_transpose(L), Tbc_t);
    relative_L <- cholesky_decompose(Lambda);
    zeros <- rep_vector(0, pm1);
    for (i in 1:N) {
      vector[pm1] u;
      int biggest;
      u <- multi_normal_cholesky_rng(zeros, relative_L);
      if (rows(alpha) > 0)     u <- u + alpha;
      if (rows(beta) > 0)      u <- u + beta * X[i];
      for (j in 1:rows(gamma)) u <- u + gamma[j] * Z[i,j];
      biggest <- sort_indices_desc(u)[1];
      if (biggest < 0) mean_PPD[1] <- mean_PPD[1] + 1;
      else mean_PPD[biggest + 1] <- mean_PPD[biggest + 1] + 1;
    }
    mean_PPD <- mean_PPD / N;
  }
}
