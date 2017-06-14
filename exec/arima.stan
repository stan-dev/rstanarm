functions {
  /*
   * Transform from partial autocorrelations to autocorrelations that
   * imply stationarity; also works with moving average terms. See
   * M.C. Jones, 1987, "Randomly Choosing Parameters from the Stationarity
   * and Invertibility Region of Autoregressive-Moving Average Models",
   * Applied Statistics, 36(2) p. 134 -- 138
   *
   * @param r vector of partial autocorrelations
   * @return vector of autocorrelations
   */
  vector LevinsonDurbin(vector r) {
    vector[rows(r)] phi;
    vector[rows(r)] work;
    phi[1] = r[1];
    work[1] = r[1];
    for (j in 2:rows(r)) {
      real a;
      a = r[j];
      if (a < -1) reject("partial autocorrelations must be >= -1");
      if (a >  1) reject("partial autocorrelations must be <= 1");
      for (k in 1:(j - 1)) work[k] = work[k] - a * phi[j - k];
      work[j] = a;
      phi[1:j] = work[1:j];
    }
    return phi;
  }
}
data {
  real lb;
  real<lower=lb> ub;
  int<lower=0> T;
  int<lower=0, upper=T> p;
  int<lower=0, upper=T> q;
  vector<lower=lb,upper=ub>[T] yy;
  matrix[T, p] X;
  // delta_{AR,MA} = 0 -> uniform prior over the stationary region
  real<lower=0> delta_AR;
  real<lower=0> delta_MA;
  int<lower=0,upper=1> has_intercept;
}
transformed data {
  vector<lower=0>[p] ARshape_1;
  vector<lower=0>[p] ARshape_2;
  vector<lower=0>[q] MAshape_1;
  vector<lower=0>[q] MAshape_2;
  for (k in 1:p) {
    real pp1mk;
    pp1mk = p + 1 - k;
    ARshape_1[k] = floor(0.5 * (k + 1)) + delta_AR * pp1mk;
    ARshape_2[k] = floor(0.5 * k) + 1   + delta_AR * pp1mk;
  }
  for (k in 1:q) {
    real qp1mk;
    qp1mk = q + 1 - k;
    MAshape_1[k] = floor(0.5 * (k + 1)) + delta_MA * qp1mk;
    MAshape_2[k] = floor(0.5 * k) + 1   + delta_MA * qp1mk;
  }
}
parameters {
  real<lower=lb, upper=ub> mu[has_intercept];
  vector<lower=-1,upper=1>[p] r_phi;
  vector<lower=-1,upper=1>[q] r_theta;
  real<lower=0> sigma;
}
transformed parameters {
  real alpha[has_intercept];
  vector[p] phi;
  vector[q] theta;
  if (p > 0) phi   = LevinsonDurbin(r_phi);
  if (q > 0) theta = LevinsonDurbin(r_theta);
  if (has_intercept) {
    if (p > 0) alpha[1] = mu[1] * (1 - sum(phi));
    else alpha = mu;
  }
}
model {
  vector[T] eta;
  if (has_intercept) {
    if (p > 0) eta = alpha[1] + X * phi;
    else       eta = rep_vector(alpha[1], T);
  }
  else {
    if (p > 0) eta = X * phi;
    else eta = rep_vector(0.0, T);
  }
  if (q > 0) {
    vector[T] epsilon;
    epsilon = yy - eta;
    for (t in 1:q) for (k in 1:t)
      eta[t] = eta[t] + theta[k] * epsilon[t - k];
    for (t in (q + 1):T) for (k in 1:q)
      eta[t] = eta[t] + theta[k] * epsilon[t - k];
  }
  target+= normal_lpdf(yy | eta, sigma);

  // ignore Jacobian warnings because derivatives are constants
  target+= beta_lpdf(0.5 * r_phi + 0.5 | ARshape_1, ARshape_2);
  target+= beta_lpdf(0.5 * r_theta + 0.5 | MAshape_1, MAshape_2);
}
