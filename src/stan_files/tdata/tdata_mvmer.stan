  // dimensions for hs priors
  int<lower=0> yHs1 = get_nvars_for_hs(M > 0 ? y_prior_dist[1] : 0);
  int<lower=0> yHs2 = get_nvars_for_hs(M > 1 ? y_prior_dist[2] : 0);
  int<lower=0> yHs3 = get_nvars_for_hs(M > 2 ? y_prior_dist[3] : 0);
  int<lower=0> yHs4 = get_nvars_for_hs(M > 3 ? y_prior_dist[4] : 0);
  int<lower=0> yHs5 = get_nvars_for_hs(M > 4 ? y_prior_dist[5] : 0);
  int<lower=0> yHs6 = get_nvars_for_hs(M > 5 ? y_prior_dist[6] : 0);
  int<lower=0> yHs7 = get_nvars_for_hs(M > 6 ? y_prior_dist[7] : 0);
  int<lower=0> yHs8 = get_nvars_for_hs(M > 7 ? y_prior_dist[8] : 0);
  int<lower=0> yHs9 = get_nvars_for_hs(M > 8 ? y_prior_dist[9] : 0);
  int<lower=0> yHs10 = get_nvars_for_hs(M > 9 ? y_prior_dist[10] : 0);
  int<lower=0> yHs11 = get_nvars_for_hs(M > 10 ? y_prior_dist[11] : 0);
  int<lower=0> yHs12 = get_nvars_for_hs(M > 11 ? y_prior_dist[12] : 0);
  int<lower=0> yHs13 = get_nvars_for_hs(M > 12 ? y_prior_dist[13] : 0);
  int<lower=0> yHs14 = get_nvars_for_hs(M > 13 ? y_prior_dist[14] : 0);
  int<lower=0> yHs15 = get_nvars_for_hs(M > 14 ? y_prior_dist[15] : 0);
  int<lower=0> yHs16 = get_nvars_for_hs(M > 15 ? y_prior_dist[16] : 0);
  int<lower=0> yHs17 = get_nvars_for_hs(M > 16 ? y_prior_dist[17] : 0);
  int<lower=0> yHs18 = get_nvars_for_hs(M > 17 ? y_prior_dist[18] : 0);
  int<lower=0> yHs19 = get_nvars_for_hs(M > 18 ? y_prior_dist[19] : 0);
  int<lower=0> yHs20 = get_nvars_for_hs(M > 19 ? y_prior_dist[20] : 0);

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
  real sum_log_y4 = M > 3 && (family[4] == 2 || family[4] == 3) ?
    sum(log(yReal4)) : not_a_number();
  real sum_log_y5 = M > 4 && (family[5] == 2 || family[5] == 3) ?
    sum(log(yReal5)) : not_a_number();
  real sum_log_y6 = M > 5 && (family[6] == 2 || family[6] == 3) ?
    sum(log(yReal6)) : not_a_number();
  real sum_log_y7 = M > 6 && (family[7] == 2 || family[7] == 3) ?
    sum(log(yReal7)) : not_a_number();
  real sum_log_y8 = M > 7 && (family[8] == 2 || family[8] == 3) ?
    sum(log(yReal8)) : not_a_number();
  real sum_log_y9 = M > 8 && (family[9] == 2 || family[9] == 3) ?
    sum(log(yReal9)) : not_a_number();
  real sum_log_y10 = M > 9 && (family[10] == 2 || family[10] == 3) ?
    sum(log(yReal10)) : not_a_number();
  real sum_log_y11 = M > 10 && (family[11] == 2 || family[11] == 3) ?
    sum(log(yReal11)) : not_a_number();
  real sum_log_y12 = M > 11 && (family[12] == 2 || family[12] == 3) ?
    sum(log(yReal12)) : not_a_number();
  real sum_log_y13 = M > 12 && (family[13] == 2 || family[13] == 3) ?
    sum(log(yReal13)) : not_a_number();
  real sum_log_y14 = M > 13 && (family[14] == 2 || family[14] == 3) ?
    sum(log(yReal14)) : not_a_number();
  real sum_log_y15 = M > 14 && (family[15] == 2 || family[15] == 3) ?
    sum(log(yReal15)) : not_a_number();
  real sum_log_y16 = M > 15 && (family[16] == 2 || family[16] == 3) ?
    sum(log(yReal16)) : not_a_number();
  real sum_log_y17 = M > 16 && (family[17] == 2 || family[17] == 3) ?
    sum(log(yReal17)) : not_a_number();
  real sum_log_y18 = M > 17 && (family[18] == 2 || family[18] == 3) ?
    sum(log(yReal18)) : not_a_number();
  real sum_log_y19 = M > 18 && (family[19] == 2 || family[19] == 3) ?
    sum(log(yReal19)) : not_a_number();
  real sum_log_y20 = M > 19 && (family[20] == 2 || family[20] == 3) ?
    sum(log(yReal20)) : not_a_number();
  vector[M > 0 && family[1] == 3 ? yNobs[1] : 0] sqrt_y1;
  vector[M > 1 && family[2] == 3 ? yNobs[2] : 0] sqrt_y2;
  vector[M > 2 && family[3] == 3 ? yNobs[3] : 0] sqrt_y3;
  vector[M > 3 && family[4] == 3 ? yNobs[4] : 0] sqrt_y4;
  vector[M > 4 && family[5] == 3 ? yNobs[5] : 0] sqrt_y5;
  vector[M > 5 && family[6] == 3 ? yNobs[6] : 0] sqrt_y6;
  vector[M > 6 && family[7] == 3 ? yNobs[7] : 0] sqrt_y7;
  vector[M > 7 && family[8] == 3 ? yNobs[8] : 0] sqrt_y8;
  vector[M > 8 && family[9] == 3 ? yNobs[9] : 0] sqrt_y9;
  vector[M > 9 && family[10] == 3 ? yNobs[10] : 0] sqrt_y10;
  vector[M > 10 && family[11] == 3 ? yNobs[11] : 0] sqrt_y11;
  vector[M > 11 && family[12] == 3 ? yNobs[12] : 0] sqrt_y12;
  vector[M > 12 && family[13] == 3 ? yNobs[13] : 0] sqrt_y13;
  vector[M > 13 && family[14] == 3 ? yNobs[14] : 0] sqrt_y14;
  vector[M > 14 && family[15] == 3 ? yNobs[15] : 0] sqrt_y15;
  vector[M > 15 && family[16] == 3 ? yNobs[16] : 0] sqrt_y16;
  vector[M > 16 && family[17] == 3 ? yNobs[17] : 0] sqrt_y17;
  vector[M > 17 && family[18] == 3 ? yNobs[18] : 0] sqrt_y18;
  vector[M > 18 && family[19] == 3 ? yNobs[19] : 0] sqrt_y19;
  vector[M > 19 && family[20] == 3 ? yNobs[20] : 0] sqrt_y20;
  vector[M > 0 && family[1] == 3 ? yNobs[1] : 0] log_y1;
  vector[M > 1 && family[2] == 3 ? yNobs[2] : 0] log_y2;
  vector[M > 2 && family[3] == 3 ? yNobs[3] : 0] log_y3;
  vector[M > 3 && family[4] == 3 ? yNobs[4] : 0] log_y4;
  vector[M > 4 && family[5] == 3 ? yNobs[5] : 0] log_y5;
  vector[M > 5 && family[6] == 3 ? yNobs[6] : 0] log_y6;
  vector[M > 6 && family[7] == 3 ? yNobs[7] : 0] log_y7;
  vector[M > 7 && family[8] == 3 ? yNobs[8] : 0] log_y8;
  vector[M > 8 && family[9] == 3 ? yNobs[9] : 0] log_y9;
  vector[M > 9 && family[10] == 3 ? yNobs[10] : 0] log_y10;
  vector[M > 10 && family[11] == 3 ? yNobs[11] : 0] log_y11;
  vector[M > 11 && family[12] == 3 ? yNobs[12] : 0] log_y12;
  vector[M > 12 && family[13] == 3 ? yNobs[13] : 0] log_y13;
  vector[M > 13 && family[14] == 3 ? yNobs[14] : 0] log_y14;
  vector[M > 14 && family[15] == 3 ? yNobs[15] : 0] log_y15;
  vector[M > 15 && family[16] == 3 ? yNobs[16] : 0] log_y16;
  vector[M > 16 && family[17] == 3 ? yNobs[17] : 0] log_y17;
  vector[M > 17 && family[18] == 3 ? yNobs[18] : 0] log_y18;
  vector[M > 18 && family[19] == 3 ? yNobs[19] : 0] log_y19;
  vector[M > 19 && family[20] == 3 ? yNobs[20] : 0] log_y20;
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
  if (M > 3 && family[4] == 3) {
    sqrt_y4 = sqrt(yReal4);
    log_y4 = log(yReal4);
  }
  if (M > 4 && family[5] == 3) {
    sqrt_y5 = sqrt(yReal5);
    log_y5 = log(yReal5);
  }
  if (M > 5 && family[6] == 3) {
    sqrt_y6 = sqrt(yReal6);
    log_y6 = log(yReal6);
  }
  if (M > 6 && family[7] == 3) {
    sqrt_y7 = sqrt(yReal7);
    log_y7 = log(yReal7);
  }
  if (M > 7 && family[8] == 3) {
    sqrt_y8 = sqrt(yReal8);
    log_y8 = log(yReal8);
  }
  if (M > 8 && family[9] == 3) {
    sqrt_y9 = sqrt(yReal9);
    log_y9 = log(yReal9);
  }
  if (M > 9 && family[10] == 3) {
    sqrt_y10 = sqrt(yReal10);
    log_y10 = log(yReal10);
  }
  if (M > 10 && family[11] == 3) {
    sqrt_y11 = sqrt(yReal11);
    log_y11 = log(yReal11);
  }
  if (M > 11 && family[12] == 3) {
    sqrt_y12 = sqrt(yReal12);
    log_y12 = log(yReal12);
  }
  if (M > 12 && family[13] == 3) {
    sqrt_y13 = sqrt(yReal13);
    log_y13 = log(yReal13);
  }
  if (M > 13 && family[14] == 3) {
    sqrt_y14 = sqrt(yReal14);
    log_y14 = log(yReal14);
  }
  if (M > 14 && family[15] == 3) {
    sqrt_y15 = sqrt(yReal15);
    log_y15 = log(yReal15);
  }
  if (M > 15 && family[16] == 3) {
    sqrt_y16 = sqrt(yReal16);
    log_y16 = log(yReal16);
  }
  if (M > 16 && family[17] == 3) {
    sqrt_y17 = sqrt(yReal17);
    log_y17 = log(yReal17);
  }
  if (M > 17 && family[18] == 3) {
    sqrt_y18 = sqrt(yReal18);
    log_y18 = log(yReal18);
  }
  if (M > 18 && family[19] == 3) {
    sqrt_y19 = sqrt(yReal19);
    log_y19 = log(yReal19);
  }
  if (M > 19 && family[20] == 3) {
    sqrt_y20 = sqrt(yReal20);
    log_y20 = log(yReal20);
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
