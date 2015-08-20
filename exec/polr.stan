# GLM for an ordinal outcome with coherent priors
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
     // "logistic", "probit", "loglog", "cloglog", "cauchit"
     real p;
     if      (link == 1) p <- inv_logit(x);
     else if (link == 2) p <- Phi(x);
     else if (link == 3) p <- gumbel_cdf(x, 0, 1);
     else if (link == 4) p <- inv_cloglog(x);
     else                p <- cauchy_cdf(x, 0, 1);
     return p;
   }
  
  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer outcome variable.
   * @param eta A vector of linear predictors
   * @param cutpoints An ordered vector of cutpoints
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_polr(int[] y, vector eta, vector cutpoints, int link) {
    vector[rows(eta)] ll;
    int N;
    int J;
    N <- rows(eta);
    J <- rows(cutpoints) + 1;
    if (link < 1 || link > 5) reject("Invalid link");
    for (n in 1:N) {
           if (y[n] == 1) ll[n] <- CDF_polr(cutpoints[1] - eta[n], link);
      else if (y[n] == J) ll[n] <- 1 - CDF_polr(cutpoints[J - 1] - eta[n], link);
      else ll[n] <- CDF_polr(cutpoints[y[n]]     - eta[n], link) - 
                    CDF_polr(cutpoints[y[n] - 1] - eta[n], link);
    }
    return log(ll);
  }

  /** 
   * The standard normal inverse CDF
   *
   * Algorithm first derived in 2003 by Peter Jon Aklam at 
   * http://home.online.no/~pjacklam/notes/invnorm/
   *
   * @param p Scalar typically on the (0,1) interval
   * @return Scalar, x, on the entire extended real line, such that Phi(x) "=" p
  */
  real inv_Phi(real p) {
    real x; real q; real r; real p_low;  real p_high; real e; real u;
    real a[6]; real b[5]; real c[6]; real d[4];
    if( p <= 0 )     return(negative_infinity());
    else if (p >= 1) return(positive_infinity());
    
    a[1] <- -3.969683028665376e+01;
    a[2] <-  2.209460984245205e+02;
    a[3] <- -2.759285104469687e+02;
    a[4] <-  1.383577518672690e+02;
    a[5] <- -3.066479806614716e+01;
    a[6] <-  2.506628277459239e+00;
  
    b[1] <- -5.447609879822406e+01;
    b[2] <-  1.615858368580409e+02;
    b[3] <- -1.556989798598866e+02;
    b[4] <-  6.680131188771972e+01;
    b[5] <- -1.328068155288572e+01;
  
    c[1] <- -7.784894002430293e-03;
    c[2] <- -3.223964580411365e-01;
    c[3] <- -2.400758277161838e+00;
    c[4] <- -2.549732539343734e+00;
    c[5] <-  4.374664141464968e+00;
    c[6] <-  2.938163982698783e+00;
  
    d[1] <-  7.784695709041462e-03;
    d[2] <-  3.224671290700398e-01;
    d[3] <-  2.445134137142996e+00;
    d[4] <-  3.754408661907416e+00;
    
    p_low  <- 0.02425;
    p_high <- 0.97575; # 1 - p_low
    if ( (p_low <= p) && (p <= p_high) ) { # central region
      q <- p - 0.5;
      r <- square(q);
      x <- (((((a[1]*r+a[2])*r+a[3])*r+a[4])*r+a[5])*r+a[6])*q /
           (((((b[1]*r+b[2])*r+b[3])*r+b[4])*r+b[5])*r+1.0);      
    }
    else if (p < p_low) { # lower region
      q <- sqrt(-2.0*log(p));
      x <- (((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /
            ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1.0);            
    }
    else { # upper region
      q <- sqrt(-2.0*log1m(p));
      x <- -(((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /
             ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1.0);
    }

    # refinement via Newton step using reciprocal PDF as derivative
    e <- Phi(x) - p; // error thusfar, usually about 1e-8
    u <- e * sqrt(2.0 * pi()) * exp(0.5 * square(x));
    x <- x - u / (1.0 + 0.5 * x * u);
    return x;
  }
  
  /**
   * Map from conditional probabilities to cutpoints
   *
   * @param probabilities A J-simplex
   * @param scale A positive number
   * @param link An integer indicating the link function
   * @return A vector of length J - 1 whose elements are in increasing order
   */
  vector make_cutpoints(vector probabilities, real scale, int link) {
    vector[rows(probabilities) - 1] cutpoints;
    real running_sum;
    // links in MASS::polr() are in a different order than binomial() 
    // "logistic", "probit", "loglog", "cloglog", "cauchit"
    if (link < 1 || link > 5) reject("invalid link");
    running_sum <- 0;
    if (link == 1) for(c in 1:(rows(cutpoints))) {
      running_sum  <- running_sum + probabilities[c];
      cutpoints[c] <- logit(running_sum);
    }
    else if (link == 2) for(c in 1:(rows(cutpoints))) {
      running_sum  <- running_sum + probabilities[c];
      cutpoints[c] <- inv_Phi(running_sum);
    }
    else if (link == 3) for(c in 1:(rows(cutpoints))) {
      running_sum  <- running_sum + probabilities[c];
      cutpoints[c] <- -log(-log(running_sum));
    }
    else if (link == 4) for(c in 1:(rows(cutpoints))) {
      running_sum  <- running_sum + probabilities[c];
      cutpoints[c] <- log(-log1m(running_sum));
    }
    else for(c in 1:(rows(cutpoints))) {
      running_sum  <- running_sum + probabilities[c];
      cutpoints[c] <- tan(pi() * (running_sum - 0.5));
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
    iter <- 0;
    ystar <- not_a_number();
    if (lower >= upper) reject("lower must be less than upper");
    if      (link == 1) while(!(ystar > lower && ystar < upper))
      ystar <- logistic_rng(eta, 1);
    else if (link == 2) while(!(ystar > lower && ystar < upper))
      ystar <- normal_rng(eta, 1);
    else if (link == 3) while(!(ystar > lower && ystar < upper))
      ystar <- gumbel_rng(eta, 1);
    else if (link == 4) while(!(ystar > lower && ystar < upper))
      ystar <- log(-log1m(uniform_rng(0,1)));
    else if (link == 5) while(!(ystar > lower && ystar < upper))
      ystar <- cauchy_rng(eta, 1);
    return ystar;
  }
}
data {
  int<lower=2> J;                # number of outcome categories, which typically is > 2
  int<lower=1> N;                # number of observations
  int<lower=0> K;                # number of predictors (excluding a constant)
  matrix[N,K]  X;                # centered predictor matrix
  vector[K] xbar;                # means of the predictors
  vector<lower=0>[K] s_X;        # standard deviations of the predictors
  int<lower=1,upper=J> y[N];     # outcome
  int<lower=0,upper=1> prior_PD; # flag indicating whether to draw from the prior predictive
  
  # link function from location to linear predictor
  int<lower=1,upper=5> link;

  # weights
  int<lower=0,upper=1> has_weights; # 0 = No, 1 = Yes
  vector[N * has_weights] weights;
  
  # offset
  int<lower=0,upper=1> has_offset;  # 0 = No, 1 = Yes
  vector[N * has_offset] offset;
  
  # prior family (zero indicates no prior!!!)
  int<lower=0,upper=1> prior_dist;  # 1 = Ben
  
  # hyperparameter values
  real<lower=0> shape;
  vector<lower=0>[J] prior_counts;
  
  int<lower=0,upper=1> do_residuals;
}
transformed data {
  real<lower=shape> shapephalf;
  real<lower=0> half_K;
  int<lower=0,upper=1> is_constant;
  matrix[K,K] middle;
  shapephalf <- shape + 0.5;
  half_K <- 0.5 * K;
  is_constant <- 1;
  for (j in 1:J) if (prior_counts[j] != 1) is_constant <- 0;
  middle <- xbar * transpose(xbar);
}
parameters {
  simplex[J] pi;
  row_vector[K] z_beta;
  cholesky_factor_corr[K] L[prior_dist == 1];
  real<lower=0,upper=1>  R2[prior_dist == 1];
}
transformed parameters {
  real Delta_y;
  vector[K] beta;
  vector[J-1] cutpoints;
  if (prior_dist == 1) {
    Delta_y <- inv(sqrt(1 - R2[1]));
    if (K > 1)
      beta <- transpose(mdivide_right_tri_low(z_beta, L[1])) *
              sqrt(R2[1] / dot_self(z_beta)) ./ s_X * Delta_y;
    else beta[1] <- sqrt(R2[1]) / s_X[1] * Delta_y;
  }
  else { // prior_dist == 0
    beta <- transpose(z_beta);
    Delta_y <- sqrt(quad_form(middle, beta) + 1);
  }
  cutpoints <- make_cutpoints(pi, Delta_y, link);
}
model {
  vector[N] eta;
  if (K > 0) eta <- X * beta;
  else eta <- rep_vector(0.0, N);
  if (has_offset == 1) eta <- eta + offset;
  if (has_weights == 0 && prior_PD == 0) { # unweighted log-likelihoods
    increment_log_prob(pw_polr(y, eta, cutpoints, link));
  }
  else if (prior_PD == 0) { # weighted log-likelihoods
    increment_log_prob(dot_product(weights, pw_polr(y, eta, cutpoints, link)));
  }
  
  if (prior_dist == 1) {
    z_beta ~ normal(0, 1);
    if (K > 1) L[1] ~ lkj_corr_cholesky(shapephalf);
    R2[1] ~ beta(half_K, shape);
    if (is_constant == 0) pi ~ dirichlet(prior_counts);
  }
  /* else prior_dist is 0 and nothing is added */
}
generated quantities {
  vector[J-1] zeta;
  vector[(J > 2) * (J - 1) + 1] mean_PPD;
  vector[N * do_residuals] residuals;
  zeta <- cutpoints + dot_product(xbar, beta);
  if (J == 2) zeta <- -zeta;
  mean_PPD <- rep_vector(0,rows(mean_PPD));
  {
    vector[N] eta;
    if (K > 0) eta <- X * beta;
    else eta <- rep_vector(0.0, N);
    if (has_offset == 1) eta <- eta + offset;
    for (n in 1:N) {
      vector[J] theta;
      int y_tilde;
      real previous;
      real ystar;
      theta[1] <- CDF_polr(cutpoints[1] - eta[n], link);
      previous <- theta[1];
      for (j in 2:(J-1)) {
        real current;
        current <- CDF_polr(cutpoints[j] - eta[n], link);
        theta[j] <- current - previous;
        previous <- current;
      }
      theta[J] <- 1 - previous;
      if (J == 2) {
        mean_PPD[1] <- mean_PPD[1] + bernoulli_rng(theta[J]);
      }
      else {
        y_tilde <- categorical_rng(theta);
        mean_PPD[y_tilde] <- mean_PPD[y_tilde] + 1;
      }
      
      if (do_residuals) {
        if (y[n] == 1)
          ystar <- draw_ystar_rng(negative_infinity(), cutpoints[1], eta[n], link);
        else if (y[n] == J)
          ystar <- draw_ystar_rng(cutpoints[J - 1], positive_infinity(), eta[n], link);
        else ystar <- draw_ystar_rng(cutpoints[y[n] - 1], cutpoints[y[n]], eta[n], link);
        residuals[n] <- ystar - eta[n];
      }
    }
    mean_PPD <- mean_PPD / N;
  }
}
