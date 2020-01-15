  vector[N] eta;  // linear predictor
  if (K > 0) {
    if (dense_X) eta = X[1] * beta;
    else eta = csr_matrix_times_vector2(N, K, w_X, v_X, u_X, beta);
  }
  else eta = rep_vector(0.0, N);
  if (has_offset == 1) eta += offset_;
  if (K_smooth) eta += S * beta_smooth;
