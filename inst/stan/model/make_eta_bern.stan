  vector[N[1]] eta0;
  vector[N[2]] eta1;
  if (K > 0) {
    if (dense_X) {
      eta0 = N[1] > 0 ? X0[1] * beta : rep_vector(0.0, 0);
      eta1 = N[2] > 0 ? X1[1] * beta : rep_vector(0.0, 0);
    }
    else {
      eta0 = csr_matrix_times_vector(N[1], K, w_X0, v_X0, u_X0, beta);
      eta1 = csr_matrix_times_vector(N[2], K, w_X1, v_X1, u_X1, beta);
    }
  }
  else {
    eta0 = rep_vector(0.0, N[1]);
    eta1 = rep_vector(0.0, N[2]);
  }
  if (has_intercept == 0 && dense_X) {
    real tmp = dot_product(xbar, beta);
    if (N[1] > 0) eta0 += tmp;
    if (N[2] > 0) eta1 += tmp;
  }
  if (has_offset == 1) {
    if (N[1] > 0) eta0 += offset0;
    if (N[2] > 0) eta1 += offset1;
  }
  if (K_smooth) {
    if (N[1] > 0) eta0 += S0 * beta_smooth;
    if (N[2] > 0) eta1 += S1 * beta_smooth;
  }
  if (special_case) for (i in 1:t) {
    if (N[1] > 0) eta0 += b[V0[i]];
    if (N[2] > 0) eta1 += b[V1[i]];
  } else if (t > 0) {
    if (N[1] > 0) eta0 += csr_matrix_times_vector(N[1], q, w0, v0, u0, b);
    if (N[2] > 0) eta1 += csr_matrix_times_vector(N[2], q, w1, v1, u1, b);
  }
