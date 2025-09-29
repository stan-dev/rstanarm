#include /include/Columbia_copyright.stan
#include /include/license.stan

// GLM for a Gaussian outcome with no link function
functions {
  /**
   * Increments the log-posterior with the logarithm of a multivariate normal
   * likelihood with a scalar standard deviation for all errors
   * Equivalent to normal_lpdf(y | intercept + Q * R * beta, sigma) but faster
   * @param theta vector of coefficients (excluding intercept), equal to R * beta
   * @param b precomputed vector of OLS coefficients (excluding intercept) in Q-space
   * @param intercept scalar (assuming columns of Q have mean zero)
   * @param ybar precomputed sample mean of the outcome
   * @param SSR positive precomputed value of the sum of squared OLS residuals
   * @param sigma positive scalar for the standard deviation of the errors
   * @param N integer equal to the number of observations
   */
  real mvn_ols_qr_lpdf(vector theta, vector b, real intercept, real ybar,
                       real SSR, real sigma, int N) {
    return -0.5 * (dot_self(theta - b) + N * square(intercept - ybar) + SSR)
           / square(sigma)
           - // 0.91... is log(sqrt(2 * pi()))
           N * (log(sigma) + 0.91893853320467267);
  }
}
data {
  int<lower=0, upper=1> has_intercept; // 0 = no, 1 = yes
  int<lower=0, upper=1> prior_dist_for_intercept; // 0 = none, 1 = normal
  real<lower=0> prior_scale_for_intercept; // 0 = by CLT
  real prior_mean_for_intercept; // expected value for alpha
  int<lower=0, upper=1> prior_dist; // 0 = uniform for R^2, 1 = Beta(K/2,eta)
  int<lower=0, upper=1> prior_PD; // 0 = no, 1 = yes to drawing from the prior
  real<lower=0> eta; // shape hyperparameter

  int<lower=1> J; // number of groups
  // the rest of these are indexed by group but should work even if J = 1
  array[J] int<lower=1> N; // number of observations
  int<lower=1, upper=min(N)> K; // number of predictors
  array[J] vector[K] xbarR_inv; // vector of means of the predictors
  array[J] real ybar; // sample mean of outcome
  real center_y; // zero or sample mean of outcome
  array[J] real<lower=0> s_Y; // standard deviation of the outcome
  array[J] vector[K] Rb; // OLS coefficients
  array[J] real<lower=0> SSR; // OLS sum-of-squared residuals
  array[J] matrix[K, K] R_inv; // inverse R matrices
}
transformed data {
  real half_K = 0.5 * K;
  array[J] real sqrt_inv_N;
  array[J] real sqrt_Nm1;
  for (j in 1 : J) {
    sqrt_inv_N[j] = sqrt(1.0 / N[j]);
    sqrt_Nm1[j] = sqrt(N[j] - 1.0);
  }
}
parameters {
  // must not call with init="0"
  // https://github.com/stan-dev/rstanarm/issues/603#issuecomment-1785928224
  array[K > 1 ? J : 0] unit_vector[K > 1 ? K : 2] u; // primitives for coefficients
  array[J * has_intercept] real z_alpha; // primitives for intercepts
  array[J] real<lower=(K > 1 ? 0 : -1), upper=1> R2; // proportions of variance explained
  vector[J * (1 - prior_PD)] log_omega; // under/overfitting factors
}
transformed parameters {
  array[J * has_intercept] real alpha; // uncentered intercepts
  array[J] vector[K] theta; // coefficients in Q-space
  array[J] real<lower=0> sigma; // error standard deviations
  for (j in 1 : J) {
    // marginal standard deviation of outcome for group j
    real Delta_y = prior_PD == 0 ? s_Y[j] * exp(log_omega[j]) : 1;

    // coefficients in Q-space
    if (K > 1)
      theta[j] = u[j] * sqrt(R2[j]) * sqrt_Nm1[j] * Delta_y;
    else
      theta[j][1] = R2[j] * sqrt_Nm1[j] * Delta_y;

    sigma[j] = Delta_y * sqrt(1 - R2[j]); // standard deviation of errors

    if (has_intercept == 1) {
      if (prior_dist_for_intercept == 0)  // no information
        alpha[j] = z_alpha[j];
      else if (prior_scale_for_intercept == 0)  // central limit theorem
        alpha[j] = z_alpha[j] * Delta_y * sqrt_inv_N[j]
                   + prior_mean_for_intercept;
      else // arbitrary informative prior

        alpha[j] = z_alpha[j] * prior_scale_for_intercept
                   + prior_mean_for_intercept;
    }
  }
}
model {
  if (prior_PD == 0)
    for (j in 1 : J) {
      // likelihood contribution for each group
      real shift = dot_product(xbarR_inv[j], theta[j]);
      target += mvn_ols_qr_lpdf(theta[j] | Rb[j], has_intercept == 1
                                                  ? alpha[j] + shift : shift, ybar[j], SSR[j], sigma[j], N[j]);
      // implicit: u[j] is uniform on the surface of a hypersphere
    }
  if (has_intercept == 1 && prior_dist_for_intercept > 0)
    target += normal_lpdf(z_alpha | 0, 1);
  if (prior_dist == 1) {
    if (K > 1)
      target += beta_lpdf(R2 | half_K, eta);
    else {
      array[J] real R2_abs = abs(R2);
      target += beta_lpdf(square(R2) | half_K, eta) + sum(log(R2_abs));
    }
  }
  // implicit: log_omega is uniform over the real line for all j
}
generated quantities {
  array[J] real mean_PPD;
  array[J] vector[K] beta;
  for (j in 1 : J) {
    real shift = dot_product(xbarR_inv[j], theta[j]);
    mean_PPD[j] = normal_rng(has_intercept == 1 ? alpha[j] + shift : shift,
                             sigma[j] * sqrt_inv_N[j]);
    beta[j] = R_inv[j] * theta[j];
  }
}
