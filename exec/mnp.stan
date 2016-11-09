functions {
  /**
  * linear predictor for the first p - 1 choices in a multinomial probit model
  *
  * @param alpha row_vector of intercepts (may be zero)
  * @param theta matrix of coefficients on centered unit-specific predictors
  * @param gamma vector of coefficients on centered alternative-specific predictors
  * @param Q matrix of centered unit-specific predictors
  * @param Z array of matrices of centered alternative-specific predictors
  * @return matrix of linear predictors
  */
  vector[] make_mu(vector alpha, vector theta1, vector[] theta,
                   matrix r, matrix Q1, matrix Q, matrix[] Z) {
    vector[rows(Q)] mu[cols(alpha)];
    if (cols(r) > 0) {
      vector[cols(r)]  gamma;
      if (cols(r) > rows(theta1)) gamma = r * tail(theta1, cols(r));
      else gamma = r * theta1;
      for (j in 1:cols(mu)) mu[j] = alpha[j] + Z[j] * gamma[j];
    }
    else for (j in 1:cols(alpha)) mu[j] = rep_vector(alpha[j], rows(Q));
    if (rows(theta) > 0) {
      mu[1] = mu[1] + Q1 * head(theta1, cols(Q1));
      for (j in 2:ncol(mu)) mu[j] = mu[j] + Q * theta[j - 1];
    }
    return mu;
  }
  
  /**
  * Create matrix of utilities for the first p - 1 choices
  *
  * The identifying assumption here is that utility sums to zero; as in
  * http://www.burgette.org/sMNP-R0-FINAL.pdf
  *
  * @param pi array of simplexes
  * @param log_best vector of maximal utility in log form
  * @param y integer array of observed choices
  * @param matrix of utilities
  */
  vector[] make_U(vector[] pi, vector log_best, int[] y) {
    matrix[size(y), rows(pi[1])] U;
    int N;
    int pm1;
    int p;

    N = size(y);
    pm1 = rows(U);
    p = pm1 + 1;
    for (i in 1:N) {
      int y_i;
      vector pi_i;
      real scale;
      real utility_best;
      y_i = y[i];
      pi_i = pi[i];
      utility_best = exp(log_best[i]);
      scale = utility_best * p;
      // this enforces the sum-to-zero constraint for p-dimensional utility
      for (j in 1:(y_i - 1)) U[i,j]   = utility_best - scale * pi_i[j];
      if (y_i < p)           U[i,y_i] = utility_best;
      for (j in (y_i+1):pm1) U[i,j]   = utility_best - scale * pi_i[j-1];
    }
    return U;
  }
  
  /** 
  * log-likelihood for a multinomial probit model
  *
  * The outcome data do not identify the location or scale of the utilities.
  * To pin down the location, we require that utility for each person sums to zero; see
  * http://www.burgette.org/sMNP-R0-FINAL.pdf
  * Thus, the utility of the last choice is excluded here to avoid singularities and
  * these utilities are multivariate normal in a subspace of dimension p - 1.
  * However, you can still multiply all the utilities by a positive constant, do
  * the same for the expectations and the zero-sum errors, and not change the
  * likelihood of the observed choices. Thus, we require that the diagonal of the
  * variance-covariance matrix of the raw errors be known, which is the same as in
  * https://scholar.google.com/scholar?cluster=11094876431829517046&hl=en&as_sdt=0,33
  *
  * @param U matrix of utilities for the first p - 1 choices
  * @param mu matrix of linear predictors for the first p - 1 choices
  * @param Lambda correlation matrix for all p unrestricted errors
  * @param Ts_t matrix that enforces the zero-sum constraint on the errors
  * @return A scalar log-likelihood
  */
  real ll_mnp_lp(vector[] U, vector[] mu, matrix Lambda, matrix Ts_t) {
    target += multi_normal(U, mu, quad_form(Lambda, Ts_t));
    return target();
  }
}
data {
  int<lower=1> N;                      // number of observations
  int<lower=0> K;                      // number of unit-specific predictors
  int<lower=0> q;                      // number of alternative-specific predictors
  int<lower=3> p;                      // number of choices
  matrix[N,K + q] Q1;                  // centered unit-specific predictors for choice 1
  matrix[N,p - 1] Z[q];                // centered alternative-specific predictors
  int<lower=1,upper=p> y[N];           // observed choice outcomes
  
  int<lower=0,upper=1> normalization;  // 0 -> topleft == 1, 1 -> trace == p - 1
  real<lower=0> regularization;
  matrix[K + q, K + q] R1_inv;
}
transformed data {
  matrix[N,K] Q;                       // centered unit-specific predictors for choices > 1
  matrix[q,q] R_inv;
  matrix[p,p] Ts;                      // transformation matrix that enforces the sum-to-zero constraint
  matrix[p,p - 1] Ts_t;                  
  matrix[p,p - 1] Tbc_t;               // transformation matrix relative to baseline category
  int pm1;
  matrix[q,q] r;
  if (K > 0) {
    Q = Q1[,1:K];
    R_inv = R1_inv[1:K,1:K];
  }
  pm1 = p - 1;
  Ts = rep_matrix(-1.0 / pm1, p, p);
  for (j in 1:pm1) Ts[j,j] = 1;
  Ts_t = Ts[,1:pm1];
  Tbc_t = append_row(rep_row_vector(-1, pm1), 
                     diag_matrix(rep_vector(1, pm1)));
  for (i in 1:q) for (j in 1:q) r[i,j] = R1_inv[K + i, K + j];
}
parameters {
  vector[pm1] alpha;                   // intercepts
  vector[K + q] theta1;
  vector[K, p - 2] theta;              // coefficients on unit-specific predictors

  cholesky_factor_corr[p] L;           // Cholesky factor of full error correlation matrix; see
  simplex[pm1] pi[N];                  // utility gap for non-best choices relative to best choice
  vector[N] log_best;
}
model {
  real dummy;
  matrix[p,p] Lambda;
  Lambda = multiply_lower_tri_self_transpose(L);
  for (i in 1:p) Lambda[i,i] = 1; // already 1 to numerical tolerace but reduces the autodiff
  
  dummy = ll_mnp_lp(make_U(pi, log_scale, y),
                    make_mu(alpha, theta1, theta, gamma, Q, Z), Lambda, Ts_t);
  // priors
  L ~ lkj_corr_cholesky(regularization);
}
generated quantities {
  vector[pm1] alpha_out;
  matrix[pm1, K] beta_out;
  vector[q] gamma_out;
  matrix[pm1,pm1] Sigma;
  vector[p] mean_PPD;                  // average posterior predictive distribution
  {
    vector[K + q] temp;
    temp = R_inv * theta1;
    if (q > 0) {
      gamma_out = tail(temp, q);
      beta_out[1,] = transpose(head(temp, K));
    }
    else if (q > 0) gamma_out = temp;
    else beta_out[1,] = temp;
    for (j in 2:pm1) beta_out[j,] = transpose(Q * theta1[j - 1]);
  }
  {
    real scale_factor;
    vector[K] beta_1;
    Sigma = quad_form(Lambda, Tbc_t);
    if (normalization == 0) scale_factor = Sigma[1,1];
    else scale_factor = pm1 / trace(Sigma);
    Sigma = Sigma / scale_factor;
    scale_factor = sqrt(scale_factor);

    alpha_out[pm1] = -sum(alpha) - alpha[1];
    beta_1 = beta_out[1,];
    for (k in 1:K) beta_out[pm1,k] = -sum(col(beta_out,k)) - beta_1[k];
    for (j in 2:pm1) {
      alpha_out[j-1] = alpha_out[j] - alpha[1];
      for (k in 1:K) beta_out[j-1,k] = beta_out[j,k] - beta_1[k];
    }
  }
  {
    matrix[p,p] TsL;
    vector[N] mu[pm1];
    TsL = Ts * L;
    mu = make_mu(alpha, theta1, theta, r, Q1, Q, Z);
    mean_PPD = rep_vector(0, p);
    for (i in 1:N) {
      vector[p] z;
      vector[p] w;
      vector[pm1] mu_i;
      vector[1] neg_sum_mu_i;
      int biggest;
      for (j in 1:p) z[j] = normal_rng(0,1);
      mu_i = mu[i];
      neg_sum_mu_i[1] = -sum(mu_i);
      w = append_row(mu_i, neg_sum_mu_i) + TsL * z;
      biggest = sort_indices_desc(w)[1];
      mean_PPD[biggest] = mean_PPD[biggest] + 1;
    }
    mean_PPD = mean_PPD / N;
  }
}
