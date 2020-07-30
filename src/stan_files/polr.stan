#include /pre/Columbia_copyright.stan
#include /pre/license.stan

// GLM for an ordinal outcome with coherent priors
functions {
  
  /** 
  * Evaluate a given CDF
  *
  * @param x The point to evaluate the CDF_polr at
  * @param link An integer indicating the link function
  * @return A scalar on (0,1)
  */
  real CDF_polr(real x, int link) {
    // links in MASS::polr() are in a different order than binomial() 
    // logistic, probit, loglog, cloglog, cauchit
    if (link == 1) return(inv_logit(x));
    else if (link == 2) return(Phi(x));
    else if (link == 3) return(gumbel_cdf(x, 0, 1));
    else if (link == 4) return(inv_cloglog(x));
    else if (link == 5) return(cauchy_cdf(x, 0, 1));
    else reject("Invalid link");
    return x; // never reached
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer outcome variable.
  * @param eta A vector of linear predictors
  * @param cutpoints An ordered vector of cutpoints
  * @param link An integer indicating the link function
  * @return A vector of log-likelihods
  */
  vector pw_polr(int[] y, vector eta, vector cutpoints, 
                 int link, real alpha) {
    int N = rows(eta);
    int J = rows(cutpoints) + 1;
    vector[N] ll;
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (alpha == 1) for (n in 1:N) {
      if (y[n] == 1) ll[n] = CDF_polr(cutpoints[1] - eta[n], link);
      else if (y[n] == J) ll[n] = 1 - CDF_polr(cutpoints[J - 1] - eta[n], link);
      else ll[n] = CDF_polr(cutpoints[y[n]]     - eta[n], link) - 
                   CDF_polr(cutpoints[y[n] - 1] - eta[n], link);
    }
    else for (n in 1:N) {
      if (y[n] == 1) ll[n] = CDF_polr(cutpoints[1] - eta[n], link) ^ alpha;
      else if (y[n] == J) ll[n] = 1 - CDF_polr(cutpoints[J - 1] - eta[n], link) ^ alpha;
      else reject("alpha not allowed with more than 2 outcome categories");
    }
    return log(ll);
  }
  
  /**
  * Map from conditional probabilities to cutpoints
  *
  * @param probabilities A J-simplex
  * @param scale A positive scalar
  * @param link An integer indicating the link function
  * @return A vector of length J - 1 whose elements are in increasing order
  */
  vector make_cutpoints(vector probabilities, real scale, int link) {
    int C = rows(probabilities) - 1; 
    vector[C] cutpoints;
    real running_sum = 0;
    // links in MASS::polr() are in a different order than binomial() 
    // logistic, probit, loglog, cloglog, cauchit
    if (link == 1) for(c in 1:C) {
      running_sum += probabilities[c];
      cutpoints[c] = logit(running_sum);
    }
    else if (link == 2) for(c in 1:C) {
      running_sum += probabilities[c];
      cutpoints[c] = inv_Phi(running_sum);
    }
    else if (link == 3) for(c in 1:C) {
      running_sum += probabilities[c];
      cutpoints[c] = -log(-log(running_sum));
    }
    else if (link == 4) for(c in 1:C) {
      running_sum += probabilities[c];
      cutpoints[c] = log(-log1m(running_sum));
    }
    else if (link == 5) for(c in 1:C) {
      running_sum += probabilities[c];
      cutpoints[c] = tan(pi() * (running_sum - 0.5));
    }
    else reject("invalid link");
    return scale * cutpoints;
  }
  
  /**
   * Randomly draw a value for utility
   *
   * @param low A scalar lower bound
   * @param high A scalar upper bound
   * @param eta A scalar linear predictor
   * @param link An integer indicating the link function
   * @return A scalar from the appropriate conditional distribution
   */
  real draw_ystar_rng(real low, real high, real eta, int link) {
    int iter = 0;
    real ystar = not_a_number();
    if (low >= high) 
      reject("low must be less than high");
    
    // links in MASS::polr() are in a different order than binomial() 
    // logistic, probit, loglog, cloglog, cauchit
    if (link == 1) while(!(ystar > low && ystar < high))
      ystar = logistic_rng(eta, 1);
    else if (link == 2) while(!(ystar > low && ystar < high))
      ystar = normal_rng(eta, 1);
    else if (link == 3) while(!(ystar > low && ystar < high))
      ystar = gumbel_rng(eta, 1);
    else if (link == 4) while(!(ystar > low && ystar < high))
      ystar = log(-log1m(uniform_rng(0,1)));
    else if (link == 5) while(!(ystar > low && ystar < high))
      ystar = cauchy_rng(eta, 1);
    else reject("invalid link");
    return ystar;
  }
  
  /** 
  * faster version of csr_matrix_times_vector
  * declared here and defined in C++
  *
  * @param m Integer number of rows
  * @param n Integer number of columns
  * @param w Vector (see reference manual)
  * @param v Integer array (see reference manual)
  * @param u Integer array (see reference manual)
  * @param b Vector that is multiplied from the left by the CSR matrix
  * @return A vector that is the product of the CSR matrix and b
  */
  vector csr_matrix_times_vector2(int m, int n, vector w, int[] v, int[] u, vector b);
}
data {
  // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x
#include /data/NKX.stan
  int<lower=2> J;             // number of outcome categories, which typically is > 2
  int<lower=1,upper=J> y[N];  // ordinal outcome
  // declares prior_PD, has_intercept, link, prior_dist, prior_dist_for_intercept
#include /data/data_glm.stan
  // declares has_weights, weights, has_offset, offset
#include /data/weights_offset.stan

  // hyperparameter values
  real<lower=0> regularization;
  vector<lower=0>[J] prior_counts;
  int<lower=0,upper=1> is_skewed;
  real<lower=0> shape;
  real<lower=0> rate;
  int<lower=0,upper=1> do_residuals;
}
transformed data {
  real<lower=0> half_K = 0.5 * K;
  real<lower=0> sqrt_Nm1 = sqrt(N - 1.0);
  int<lower=0,upper=1> is_constant = 1;
  vector[0] beta_smooth;  // not used
  for (j in 1:J) if (prior_counts[j] != 1) is_constant = 0;
}
parameters {
  simplex[J] pi;
  unit_vector[K] u[K > 1];
  real<lower=(K > 1 ? 0 : -1),upper=1> R2;
  real<lower=0> alpha[is_skewed];
}
transformed parameters {
  vector[K] beta;
  vector[J-1] cutpoints;
  {
    real Delta_y; 
    if (K > 1) {
      Delta_y = inv_sqrt(1 - R2);
      beta = u[1] * sqrt(R2) * Delta_y * sqrt_Nm1;
    }
    else {
      Delta_y = inv_sqrt(1 - square(R2));
      beta[1] = R2 * Delta_y * sqrt_Nm1;
    }
    cutpoints = make_cutpoints(pi, Delta_y, link);
  }
}
model {
  if (prior_PD == 0) {
#include /model/make_eta.stan
    if (has_weights == 0) {  // unweighted log-likelihoods
      if (is_skewed == 0)
        target += pw_polr(y, eta, cutpoints, link, 1.0);
      else target += pw_polr(y, eta, cutpoints, link, alpha[1]);
    }
    else {  // weighted log-likelihoods
      if (is_skewed == 0)
        target += dot_product(weights, pw_polr(y, eta, cutpoints, link, 1.0));
      else target += dot_product(weights, pw_polr(y, eta, cutpoints, link, alpha[1]));
    }
  }

  if (is_constant == 0) target += dirichlet_lpdf(pi | prior_counts);
  // implicit: u is uniform on the surface of a hypersphere
  if (prior_dist == 1) {
    if (K > 1) target += beta_lpdf( R2  | half_K, regularization);
    else target += beta_lpdf(square(R2) | half_K, regularization) + log(fabs(R2));
  }
  if (is_skewed == 1)  target += gamma_lpdf(alpha | shape, rate);
}
generated quantities {
  vector[J > 2 ? J : 1] mean_PPD = rep_vector(0, J > 2 ? J : 1);
  vector[do_residuals ? N : 0] residuals;
  vector[J-1] zeta;
  
  // xbar is actually post multiplied by R^-1
  if (dense_X) zeta = cutpoints + dot_product(xbar, beta);
  else zeta = cutpoints;
  if (J == 2) zeta *= -1.0;
  {
#include /model/make_eta.stan
    for (n in 1:N) {
      int y_tilde;
      vector[J] theta;
      real previous;
      real first = CDF_polr(cutpoints[1] - eta[n], link);
      previous = first;
      if (is_skewed) theta[1] = first ^ alpha[1];
      else theta[1] = first;
      for (j in 2:(J-1)) {
        real current = CDF_polr(cutpoints[j] - eta[n], link);
        theta[j] = current - previous;
        previous = current;
      }
      if (is_skewed == 0) theta[J] = 1 - previous;
      else theta[J] = 1 - previous ^ alpha[1];
      if (previous <= 0 || previous >= 1) {
        // do nothing
      }
      else if (J == 2) {
        mean_PPD[1] += bernoulli_rng(theta[J]);
      }
      else {
        y_tilde = categorical_rng(theta);
        mean_PPD[y_tilde] += 1;
      }
      
      if (do_residuals) {
        real ystar;
        if (y[n] == 1)
          ystar = draw_ystar_rng(negative_infinity(), cutpoints[1], eta[n], link);
        else if (y[n] == J)
          ystar = draw_ystar_rng(cutpoints[J - 1], positive_infinity(), eta[n], link);
        else ystar = draw_ystar_rng(cutpoints[y[n] - 1], cutpoints[y[n]], eta[n], link);
        residuals[n] = ystar - eta[n];
      }
    }
    mean_PPD /= (N + 0.0);
  }
}
