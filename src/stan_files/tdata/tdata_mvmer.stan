  // dimensions for hs priors
  int<lower=0> yHs1 = get_nvars_for_hs(M > 0 ? y_prior_dist[1] : 0);
  int<lower=0> yHs2 = get_nvars_for_hs(M > 1 ? y_prior_dist[2] : 0);
  int<lower=0> yHs3 = get_nvars_for_hs(M > 2 ? y_prior_dist[3] : 0);

  // data for decov prior
  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p) * (t > 0);
  int<lower=0> len_rho = sum(p) - t;
  int<lower=1> pos = 1;
  array[len_concentration] real<lower=0> delta;

  // data for lkj prior
  array[prior_dist_for_cov == 2 ? (bK1 + choose(bK1, 2)) : 0] int bCov1_idx;
  array[prior_dist_for_cov == 2 ? (bK2 + choose(bK2, 2)) : 0] int bCov2_idx;

  // transformations of data
  real sum_log_y1 = M > 0 && (family[1] == 2 || family[1] == 3) ?
    sum(log(yReal1)) : not_a_number();
  real sum_log_y2 = M > 1 && (family[2] == 2 || family[2] == 3) ?
    sum(log(yReal2)) : not_a_number();
  real sum_log_y3 = M > 2 && (family[3] == 2 || family[3] == 3) ?
    sum(log(yReal3)) : not_a_number();
  vector[M > 0 && family[1] == 3 ? yNobs[1] : 0] sqrt_y1;
  vector[M > 1 && family[2] == 3 ? yNobs[2] : 0] sqrt_y2;
  vector[M > 2 && family[3] == 3 ? yNobs[3] : 0] sqrt_y3;
  vector[M > 0 && family[1] == 3 ? yNobs[1] : 0] log_y1;
  vector[M > 1 && family[2] == 3 ? yNobs[2] : 0] log_y2;
  vector[M > 2 && family[3] == 3 ? yNobs[3] : 0] log_y3;
  if (M > 0 && family[1] == 3) {
    sqrt_y1 = sqrt(yReal1);
    log_y1 = log(yReal1);
  }
  if (M > 1 && family[2] == 3) {
    sqrt_y2 = sqrt(yReal2);
    log_y2 = log(yReal2);
  }
  if (M > 2 && family[3] == 3) {
    sqrt_y3 = sqrt(yReal3);
    log_y3 = log(yReal3);
  }

  // data for decov prior
  if (prior_dist_for_cov == 1) {
    for (i in 1:t) {
      if (p[i] > 1) {
        for (j in 1:p[i]) {
          delta[pos] = b_prior_concentration[j];
          pos += 1;
        }
      }
      for (j in 3:p[i]) len_z_T += p[i] - 1;
    }
  }

  // data for lkj prior
  if (prior_dist_for_cov == 2) {
    if (bK1 > 0)
      bCov1_idx = lower_tri_indices(bK1);
    if (bK2 > 0)
      bCov2_idx = lower_tri_indices(bK2);
  }
