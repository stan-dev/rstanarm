#include "license.stan"

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
  real ll_mvn_ols_qr_lp(vector theta, vector b,
                        real intercept, real ybar,
                        real SSR, real sigma, int N) {
    target += -0.5 * (dot_self(theta - b) + 
      N * square(intercept - ybar) + SSR) / 
      square(sigma) -// 0.91... is log(sqrt(2 * pi()))
      N * (log(sigma) + 0.91893853320467267);
    return target();
  }

  /**
   * Approximation of lgamma(z); see
   *
   * https://en.wikipedia.org/wiki/Stirling%27s_approximation#Versions_suitable_for_calculators
   *
   * @param z Positive scalar
   */
  real Nemes(real z)
    return 0.5 * log(6.283185307179586232 / z) + z * (log(z + 1 / (12 * z - 0.1 / z)) - 1);
      
  int signum(real x)
    return x < 0 ? -1 : x > 0;

  /**
   * Customized version of the bisection method to find when log-sum-exp can terminate; see
   *
   * https://en.wikipedia.org/wiki/Bisection_method#Algorithm
   *
   * This is only called from the log_besselI() function below
   *
   * @param v Order for log_besselI() function
   * @param biggest Pivot in the log-sum-exp
   * @param biggest_m Integer indicating what iteration reaches the pivot
   * @param log_half_x Real number equal to log(0.5 * x)
   */
  int bisection(real v, real biggest, int biggest_m, real log_half_x) {
    int c;
    int a;
    int b;
    int N;
    real smallest;
    real val_a;
    a = biggest_m;
    b = 2147483647; // maximum integer
    N = 1;
    smallest = -745.13321910194122211;
    val_a = -smallest;
    while (N <= 100) {
      real val;
      c = a / 2 + b / 2;
      if (a == (b - 1)) return c;
      else {
        val = -2 * Nemes(c) + (2 * c + v) * log_half_x - biggest - smallest;
        if (round(val) == 0) return c;
        else if (signum(val) == signum(val_a)) {
          a = c;
          val_a = val;
        }
        else b = c;
        N = N + 1;
      }
    }
    return c;
  }
  
  /**
   * Modded Bessel function of the first kind in log units
   *
   * This implementation is designed to not overflow or eat too much RAM
   *
   * @param x Real non-negative scalar
   * @param v Real non-negative order but may be integer
   */
  real log_besselI(real x, real v) {
    real log_half_x;
    real lfac_0;
    real lfac;
    real lgam;
    real lgam_0;
    real lcons;
    real lcons_0;
    real biggest;
    real piece;
    real smallest;
    int m;
    int cutoff;
    if (x < 0) reject("x is assumed to be non-negative")
    if (v < 0) reject("v is assumed to be non-negative");
    if (x > 700) return x - 0.5 * log(2 * pi() * x);
    log_half_x = log(0.5 * x);
    lfac_0 = 0;
    lfac = lfac_0;
    lgam_0 = lgamma(v + 1);
    lgam = lgam_0;
    lcons_0 = v * log_half_x;
    lcons = lcons_0;
    biggest = -lgam + lcons;
    piece = biggest;
    smallest = -745.13321910194122211; // exp(smallest) > 0 minimally
    m = 1;
    while (piece >= biggest) { // find pivot for log-sum-exp
      biggest = piece;
      lfac = lfac + log(m);
      lgam = lgam + log(v + m);
      lcons = lcons + 2 * log_half_x;
      piece = -lfac - lgam + lcons;
      m = m + 1;
    }
    // find next m such that exp(...) = 0
    cutoff = bisection(v, biggest, m - 1, log_half_x);
    { // start over from m = 0 and fill up summands
      vector[cutoff + 1] summands;
      lfac = lfac_0;
      lgam = lgam_0;
      lcons = lcons_0;
      piece = -lgam + lcons - biggest;
      summands[1] = exp(piece);
      m = 1;
      while (m <= cutoff) {
        lfac = lfac + log(m);
        lgam = lgam + log(v + m);
        lcons = lcons + 2 * log_half_x;
        piece = -lfac - lgam + lcons - biggest;
        m = m + 1;
        summands[m] = exp(piece);
      }
      return biggest + log(sum(summands));
    }
    reject ("should never reach here");
    return negative_infinity();
  }
  
  /**
   * von Mises-Fisher distribution in log units 
   * @param u J-array of unit vectors
   * @param mu Unit vector that is the expectation for u
   * @param kappa Non-negative concentration parameter
   */
  real VMF(vector[] u, vector mu, real kappa) {
    int p;
    real half_p;
    real half_pm1;
    int J;
    real out;
    p = rows(mu);
    if (p < 2) reject("p must be >= 2");
    half_p = 0.5 * p;
    J = size(u);
    if (kappa == 0) return -J * half_p * 1.8378770664093453391; 
    half_pm1 = half_p - 1;
    out = 0;
    for (j in 1:J) out = out + dot_product(mu, u[j]);
    out = out * kappa;
    if (p == 2) 
      return J * (-1.8378770664093453391 - log(modified_bessel_first_kind(0, kappa))) + out;
    if (p == 3) 
      return J * (log(kappa) - log(12.566370614359172464 * sinh(kappa))) + out;
    return J * ( half_pm1 * log(kappa) - half_p * 1.8378770664093453391 
                  - log_besselI(kappa, half_pm1) ) + out;
  }
}
data {
  int<lower=0,upper=1> has_intercept; // 0 = no, 1 = yes
  int<lower=0,upper=1> prior_dist_for_intercept; // 0 = none, 1 = normal
  real<lower=0> prior_scale_for_intercept;       // 0 = by CLT
  real prior_mean_for_intercept;      // expected value for alpha
  int<lower=0,upper=1> prior_dist;    // 0 = uniform for R^2, 1 = Beta(K/2,eta)
  real<lower=0> kappa_mean;           // prior expectation of concentration parameter
  int<lower=0,upper=1> prior_PD;      // 0 = no, 1 = yes to drawing from the prior
  real<lower=0> eta;                  // shape hyperparameter
  
  int<lower=1> J;                     // number of groups
  // the rest of these are indexed by group but should work even if J = 1
  int<lower=1> N[J];                  // number of observations
  int<lower=1,upper=min(N)> K;        // number of predictors
  vector[K] xbarR_inv[J];             // vector of means of the predictors
  real ybar[J];                       // sample mean of outcome
  real center_y;                      // zero or sample mean of outcome
  real<lower=0> s_Y[J];               // standard deviation of the outcome
  vector[K] Rb[J];                    // OLS coefficients
  real<lower=0> SSR[J];               // OLS sum-of-squared residuals
  matrix[K,K] R_inv[J];               // inverse R matrices
}
transformed data {
  real half_K;
  real sqrt_inv_N[J];
  real sqrt_Nm1[J];
  half_K = 0.5 * K;
  for (j in 1:J) {
    sqrt_inv_N[j] = sqrt(1.0 / N[j]);
    sqrt_Nm1[j] = sqrt(N[j] - 1.0);
  }
}
parameters { // must not call with init="0"
  unit_vector[K] u[J];                  // primitives for coefficients
  real z_alpha[J * has_intercept];      // primitives for intercepts
  real<lower=0,upper=1> R2[J];          // proportions of variance explained
  real<lower=0,upper=1> SMC[J > 1];     // proportion  of variance explained overall
  unit_vector[K] mu[J > 1];             // global location parameter
  real<lower=0> kappa[J > 1];           // unscaled concentration parameter
  vector[J * (1 - prior_PD)] log_omega; // under/overfitting factors
}
transformed parameters {
  real alpha[J * has_intercept];   // uncentered intercepts
  vector[K] theta[J];              // coefficients in Q-space
  real<lower=0> sigma[J];          // error standard deviations
  for (j in 1:J) {
    real Delta_y; // marginal standard deviation of outcome for group j
    if (prior_PD == 0) Delta_y = s_Y[j] * exp(log_omega[j]);
    else Delta_y = 1;
    
    // coefficients in Q-space
    if (K > 1) theta[j] = u[j] * sqrt(R2[j]) * sqrt_Nm1[j] * Delta_y;
    else theta[j][1] = u[j][1] * sqrt(R2[j]) * sqrt_Nm1[j] * Delta_y;
    
    sigma[j] = Delta_y * sqrt(1 - R2[j]); // standard deviation of errors
    
    if (has_intercept == 1) {
      if (prior_dist_for_intercept == 0)       // no information
        alpha[j] = z_alpha[j];
      else if (prior_scale_for_intercept == 0) // central limit theorem
        alpha[j] = z_alpha[j] * Delta_y * sqrt_inv_N[j] + prior_mean_for_intercept;
      else                                     // arbitrary informative prior
         alpha[j] = z_alpha[j] * prior_scale_for_intercept + 
                     prior_mean_for_intercept;
    }
  }
}
model {
  for (j in 1:J) { // likelihood contribution for each group
    if (prior_PD == 0) {
      real dummy; // irrelevant but useful for testing user-defined function
      real shift;
      shift = dot_product(xbarR_inv[j], theta[j]);
      dummy = ll_mvn_ols_qr_lp(theta[j], Rb[j], 
                               has_intercept == 1 ? alpha[j] + shift : shift,
                               ybar[j], SSR[j], sigma[j], N[j]);
    }
  }
  if (has_intercept == 1 && prior_dist_for_intercept > 0) 
    target += normal_lpdf(z_alpha | 0, 1);
  if (J == 1) {
    if (prior_dist == 1) target += beta_lpdf(R2 | half_K, eta);
    // implicit: u[1] is uniform on the surface of a hypersphere
  }
  else {
    target += beta_lpdf(R2 | half_K, half_K * (inv(SMC[1]) - 1));
    if (prior_dist == 1) target += beta_lpdf(SMC[1] | half_K, eta);
    target += VMF(u, mu[1], kappa_mean * kappa[1]);
    // implicit: mu is uniform on the surface of a hyperspere
    target += exponential_lpdf(kappa | 1);
  }
  // implicit: log_omega is uniform over the real line for all j
  
}
generated quantities {
  real mean_PPD[J];
  vector[K] beta[J];
  for (j in 1:J) {
    real shift;
    shift = dot_product(xbarR_inv[j], theta[j]);
    mean_PPD[j] = normal_rng(has_intercept == 1 ? alpha[j] + shift : shift,
                             sigma[j] * sqrt_inv_N[j]);
    beta[j] = R_inv[j] * theta[j];
  }
}
