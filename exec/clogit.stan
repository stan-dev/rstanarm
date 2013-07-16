# GLM for a binary outcome model --- given 1 success per group --- via Stan
data {
  int<lower=1> N; # number of observations
  int<lower=1> K; # number of predictors
  matrix[N,K]  X; # predictor matrix with no constant
  vector<lower=0,upper=1>[N] y; # outcome (usually binary, sometimes shares)
  int<lower=1> J; # number of groups
  int<lower=1> obs_per_group[J]; # number of consecutive observations per group

  int<lower=0,upper=1> has_weights; # 0 = No (weights is a ones vector), 1 = Yes
  vector[N] weights;

  int<lower=0,upper=1> has_offset;  # 0 = No (offset is a zero vector), 1 = Yes
  vector[N] offset;

  # values for hyperparameters
  vector<lower=0>[K] prior_scale;
  vector[K] prior_mean;
  vector<lower=0>[K] prior_df;
}
transformed data {
  vector[N] weighted_y;
  weighted_y <- weights .* y;
}
parameters {
  vector[K] beta;
}
model {
  int counter;
  vector[N] eta;
  counter <- 1;
  eta <- X * beta;
  if(has_offset == 1) eta <- eta + offset;

  lp__ <- lp__ + dot_product(weighted_y,eta);
  for (j in 1:J) {
    real denom[obs_per_group[j]];
    vector[obs_per_group[j]] w;
    for (k in 1:obs_per_group[j]) {
      denom[k] <- eta[counter];
      w[k] <- weighted_y[counter];
      counter <- counter + 1;
    }
    lp__ <- lp__ - sum(w * log_sum_exp(denom));
  }

  # log-priors
  for(k in 1:K) beta[k] ~ student_t(prior_df[k], prior_mean[k], prior_scale[k]);
}
