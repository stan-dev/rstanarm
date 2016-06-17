#include "license.stan"

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
    real p;
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1) p = inv_logit(x);
    else if (link == 2) p = Phi(x);
    else if (link == 3) p = gumbel_cdf(x, 0, 1);
    else if (link == 4) p = inv_cloglog(x);
    else p = cauchy_cdf(x, 0, 1);
    return p;
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
    vector[rows(eta)] ll;
    int N;
    int J;
    N = rows(eta);
    J = rows(cutpoints) + 1;
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
      else reject("alpha not allowed with more than 2 outcome categories")
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
    vector[rows(probabilities) - 1] cutpoints;
    real running_sum;
    // links in MASS::polr() are in a different order than binomial() 
    // logistic, probit, loglog, cloglog, cauchit
    if (link < 1 || link > 5) 
      reject("invalid link");
      
    running_sum = 0;
    if (link == 1) for(c in 1:(rows(cutpoints))) {
      running_sum  = running_sum + probabilities[c];
      cutpoints[c] = logit(running_sum);
    }
    else if (link == 2) for(c in 1:(rows(cutpoints))) {
      running_sum  = running_sum + probabilities[c];
      cutpoints[c] = inv_Phi(running_sum);
    }
    else if (link == 3) for(c in 1:(rows(cutpoints))) {
      running_sum  = running_sum + probabilities[c];
      cutpoints[c] = -log(-log(running_sum));
    }
    else if (link == 4) for(c in 1:(rows(cutpoints))) {
      running_sum  = running_sum + probabilities[c];
      cutpoints[c] = log(-log1m(running_sum));
    }
    else for(c in 1:(rows(cutpoints))) {
      running_sum  = running_sum + probabilities[c];
      cutpoints[c] = tan(pi() * (running_sum - 0.5));
    }
    return scale * cutpoints;
  }
  
  /**
   * Randomly draw a value for utility
   *
   * @param lower A scalar lower bound
   * @param upper A scalar upper bound
   * @param eta A scalar linear predictor
   * @param link An integer indicating the link function
   * @return A scalar from the appropriate conditional distribution
   */
  real draw_ystar_rng(real lower, real upper, real eta, int link) {
    int iter;
    real ystar;
    iter = 0;
    ystar = not_a_number();
    if (lower >= upper) 
      reject("lower must be less than upper");
    
    // links in MASS::polr() are in a different order than binomial() 
    // logistic, probit, loglog, cloglog, cauchit
    if (link < 1 || link > 5) 
      reject("invalid link");
      
    if (link == 1) while(!(ystar > lower && ystar < upper))
      ystar = logistic_rng(eta, 1);
    else if (link == 2) while(!(ystar > lower && ystar < upper))
      ystar = normal_rng(eta, 1);
    else if (link == 3) while(!(ystar > lower && ystar < upper))
      ystar = gumbel_rng(eta, 1);
    else if (link == 4) while(!(ystar > lower && ystar < upper))
      ystar = log(-log1m(uniform_rng(0,1)));
    else if (link == 5) while(!(ystar > lower && ystar < upper))
      ystar = cauchy_rng(eta, 1);
    return ystar;
  }
}
data {
  #include "NKX.stan"
  int<lower=2> J;  // number of outcome categories, which typically is > 2
  int<lower=1,upper=J> y[N];  // ordinal outcome
  #include "data_glm.stan"
  #include "weights_offset.stan"

  # hyperparameter values
  real<lower=0> regularization;
  vector<lower=0>[J] prior_counts;
  int<lower=0,upper=1> is_skewed;
  real<lower=0> shape;
  real<lower=0> rate;
  int<lower=0,upper=1> do_residuals;
}
transformed data {
  real<lower=0> half_K;
  int<lower=0,upper=1> is_constant;
  real<lower=0> sqrt_Nm1;
  half_K = 0.5 * K;
  is_constant = 1;
  for (j in 1:J) if (prior_counts[j] != 1) is_constant = 0;
  sqrt_Nm1 = sqrt(N - 1.0);
}
parameters {
  simplex[J] pi;
  vector[K] z_beta;
  real<lower=0,upper=1> R2;
  real<lower=0> alpha[is_skewed];
}
transformed parameters {
  vector[K] beta;
  vector[J-1] cutpoints;
  {
    real Delta_y;
    Delta_y = inv(sqrt(1 - R2));
    if (K > 1)
      beta = z_beta * sqrt(R2 / dot_self(z_beta)) * Delta_y * sqrt_Nm1;
    else beta[1] = sqrt(R2) * Delta_y * sqrt_Nm1;
    cutpoints = make_cutpoints(pi, Delta_y, link);
  }
}
model {
  #include "make_eta.stan"
  if (has_weights == 0 && prior_PD == 0) {  // unweighted log-likelihoods
    if (is_skewed == 0)
      target += pw_polr(y, eta, cutpoints, link, 1.0);
    else target += pw_polr(y, eta, cutpoints, link, alpha[1]);
  }
  else if (prior_PD == 0) {  // weighted log-likelihoods
    if (is_skewed == 0)
      target += dot_product(weights, pw_polr(y, eta, cutpoints, link, 1.0));
    else target += dot_product(weights, pw_polr(y, eta, cutpoints, link, alpha[1]));
  }

  if (is_constant == 0) pi ~ dirichlet(prior_counts);
  z_beta ~ normal(0,1);
  if (prior_dist == 1) R2 ~ beta(half_K, regularization);
  if (is_skewed == 1) alpha ~ gamma(shape, rate);
}
generated quantities {
  vector[J-1] zeta;
  vector[(J > 2) * (J - 1) + 1] mean_PPD;
  vector[N * do_residuals] residuals;
  
  // xbar is actually post multiplied by R^-1
  zeta = cutpoints + dot_product(xbar, beta);
  if (J == 2) zeta = -zeta;
  mean_PPD = rep_vector(0,rows(mean_PPD));
  {
    #include "make_eta.stan"
    for (n in 1:N) {
      vector[J] theta;
      int y_tilde;
      real previous;
      real ystar;
      theta[1] = CDF_polr(cutpoints[1] - eta[n], link);
      if (is_skewed) theta[1] = theta[1] ^ alpha[1];
      previous = theta[1];
      for (j in 2:(J-1)) {
        real current;
        current = CDF_polr(cutpoints[j] - eta[n], link);
        theta[j] = current - previous;
        previous = current;
      }
      if (is_skewed == 0) theta[J] = 1 - previous;
      else theta[J] = 1 - previous ^ alpha[1];
      if (previous <= 0 || previous >= 1) {
        // do nothing
      }
      else if (J == 2) {
        mean_PPD[1] = mean_PPD[1] + bernoulli_rng(theta[J]);
      }
      else {
        y_tilde = categorical_rng(theta);
        mean_PPD[y_tilde] = mean_PPD[y_tilde] + 1;
      }
      
      if (do_residuals) {
        if (y[n] == 1)
          ystar = draw_ystar_rng(negative_infinity(), cutpoints[1], eta[n], link);
        else if (y[n] == J)
          ystar = draw_ystar_rng(cutpoints[J - 1], positive_infinity(), eta[n], link);
        else ystar = draw_ystar_rng(cutpoints[y[n] - 1], cutpoints[y[n]], eta[n], link);
        residuals[n] = ystar - eta[n];
      }
    }
    mean_PPD = mean_PPD / N;
  }
}
