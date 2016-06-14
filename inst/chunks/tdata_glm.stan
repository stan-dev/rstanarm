  int<lower=0> hs;
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;
  int<lower=0,upper=1> t_any_124;
  int<lower=0,upper=1> t_all_124;
  if (prior_dist <= 2) hs <- 0;
  else if (prior_dist == 3) hs <- 2;
  else if (prior_dist == 4) hs <- 4;
  len_z_T <- 0;
  len_var_group <- sum(p) * (t > 0);
  len_rho <- sum(p) - t;
  pos <- 1;
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] <- concentration[j];
        pos <- pos + 1;
      }
    }
    for (j in 3:p[i]) len_z_T <- len_z_T + p[i] - 1;
  }
  if (prior_dist == 2) {
    t_any_124 <- 0;
    t_all_124 <- 1;
    for (k in 1:K) {
      if (prior_df[k] == 1 || prior_df[k] == 2 || prior_df[k] == 4)
        t_any_124 <- 1;
      else t_all_124 <- 0;
    }
  }
  else {
    t_any_124 <- 0;
    t_all_124 <- 0;
  }
