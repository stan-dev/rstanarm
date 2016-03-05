functions {
  /**
  * linear predictor for a multinomial probit model
  *
  * @param alpha row_vector of intercepts (may be zero)
  * @param theta matrix of coefficients on centered unit-specific predictors
  * @param gamma vector of coefficients on centered alternative-specific predictors
  * @param Q matrix of centered unit-specific predictors
  * @param Z array of matrices of centered alternative-specific predictors
  * @return matrix of linear predictors
  */
  matrix make_mu(row_vector alpha, matrix theta, vector gamma,
                 matrix Q, matrix[] Z) {
    matrix[rows(Q), cols(alpha)] mu;
    mu <- rep_matrix(alpha, rows(Q));
    for (j in 1:rows(gamma)) mu <- mu + Z[j] * gamma[j];
    if (rows(theta) > 0)     mu <- mu + Q * theta;
    return transpose(mu); // FIXME
  }
  
  /**
  * Create matrix of utilities for the first p - 1 choices
  *
  * The identifying assumption here is that utility sums to zero; as in
  * http://www.burgette.org/sMNP-R0-FINAL.pdf
  *
  * @param G array of vectors of primitives that imply zero-sum utility
  * @param y integer array of observed choices
  * @param matrix of utilities
  */
  matrix make_U(vector[] G, int[] y) {
    matrix[rows(G[1]), size(y)] U;
    int N;
    int pm1;
    int p;
    int neg_p;
    
    N <- size(y);
    pm1 <- rows(U);
    p <- pm1 + 1;
    neg_p <- -p;
    for (i in 1:N) {
      int y_i;
      vector[pm1] G_i;
      real utility_best;
      y_i <- y[i];
      G_i <- G[i];
      utility_best <- sum(G_i) / neg_p; // G[i] < 0 so utility_best > 0
      // this enforces the sum-to-zero constraint for p-dimensional utility
      for (j in 1:(y_i - 1)) U[j,i]   <- utility_best + G_i[j];
      if (y_i < p)           U[y_i,i] <- utility_best;
      for (j in (y_i+1):pm1) U[j,i]   <- utility_best + G_i[j-1];
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
  * @param Ts symmetric matrix that enforces the zero-sum constraint on the errors
  * @return A scalar log-likelihood
  */
  real ll_mnp_lp(matrix U, matrix mu, matrix Lambda, matrix Ts_t) {
    matrix[rows(U),rows(U)] L;
    L <- cholesky_decompose(quad_form(Lambda, Ts_t));
    increment_log_prob(-0.5 * sum(columns_dot_self(mdivide_left_tri_low(L, U - mu))));
    increment_log_prob(-cols(U) * sum(log(diagonal(L))));
    return get_lp();
  }
}
data {
  int<lower=1> N;                      // number of observations
  int<lower=0> K;                      // number of unit-specific predictors
  matrix[N,K] Q;                       // centered  unit-specific predictors
  int<lower=3> p;                      // number of choices
  int<lower=1,upper=p> y[N];           // observed choice outcomes
  int<lower=0> q;                      // number of alternative-specific predictors
  matrix[N,p-1] Z[q];                  // centered alternative-specific predictors
  
  int<lower=0,upper=1> normalization;  // 0 -> topleft = 1, 1 -> trace = p - 1
  real<lower=0> eta;
  matrix[p-1,p-1] R_inv;
}
transformed data {
  matrix[p,p] Ts;                      // transformation matrix that enforces the sum-to-zero constraint
  matrix[p,p-1] Ts_t;                  
  matrix[p,p-1] Tbc_t;                 // transformation matrix relative to baseline category
  int pm1;
  pm1 <- p-1;
  Ts <- rep_matrix(-1.0 / pm1, p, p);
  for (j in 1:pm1) Ts[j,j] <- 1;
  Ts_t <- Ts[,1:pm1];
  Tbc_t <- append_row(rep_row_vector(-1, pm1), 
                      diag_matrix(rep_vector(1, pm1)));  
}
parameters {
  row_vector[pm1] alpha;               // intercepts
  matrix[pm1,K] theta;                 // coefficients on unit-specific predictors
  vector[q] gamma;                     // coefficients on alternative-specific predictors
  
  cholesky_factor_corr[p] L;           // Cholesky factor of full error correlation matrix; see
  vector<upper=0>[pm1] G[N];           // utility gap for non-best choices relative to best choice
}
transformed parameters {
  matrix[p,p] Lambda;
  Lambda <- multiply_lower_tri_self_transpose(L);
  for (i in 1:p) Lambda[i,i] <- 1; // already 1 to numerical tolerace but reduces the autodiff
}
model {
  real dummy;
  dummy <- ll_mnp_lp(make_U(G, y), make_mu(alpha, theta, gamma, Q, Z), Lambda, Ts_t);
  // priors
  L ~ lkj_corr_cholesky(eta);
}
generated quantities {
  row_vector[pm1] alpha_out;
  matrix[pm1, K] beta_out;
  matrix[pm1,pm1] Sigma;
  vector[p] mean_PPD;                  // average posterior predictive distribution
  beta_out <- R_inv * theta;
  {
    real scale_factor;
    vector[K] beta_1;
    Sigma <- quad_form(Lambda, Tbc_t);
    if (normalization == 0) scale_factor <- Sigma[1,1];
    else scale_factor <- pm1 / trace(Sigma);
    Sigma <- Sigma / scale_factor;
    scale_factor <- sqrt(scale_factor);

    alpha_out[pm1] <- -sum(alpha) - alpha[1];
    beta_1 <- beta_out[,1];
    for (k in 1:K) beta_out[pm1,k] <- -sum(col(beta_out,k)) - beta_1[k];
    for (j in 2:pm1) {
      alpha_out[j-1] <- alpha_out[j] - alpha[1];
      for (k in 1:K) beta_out[k,j-1] <- beta_out[k,j] - beta_1[k];
    }
  }
  mean_PPD <- rep_vector(0, p);
  {
    matrix[p,p] TsL;
    matrix[pm1,N] mu;
    TsL <- Ts * L;
    mu <- make_mu(alpha, theta, gamma, Q, Z);
    for (i in 1:N) {
      vector[p] z;
      vector[p] w;
      vector[pm1] mu_i;
      vector[1] neg_sum_mu_i;
      int biggest;
      for (j in 1:p) z[j] <- normal_rng(0,1);
      mu_i <- col(mu, i);
      neg_sum_mu_i[1] <- -sum(mu_i);
      w <- append_row(mu_i, neg_sum_mu_i) + TsL * z;
      biggest <- sort_indices_desc(w)[1];
      mean_PPD[biggest] <- mean_PPD[biggest] + 1;
    }
    mean_PPD <- mean_PPD / N;
  }
}
