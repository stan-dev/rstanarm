  vector[N[1]] eta0;
  vector[N[2]] eta1;
  if (K > 0) {
    if (dense_X) {
      eta0 = X0[1] * beta;
      eta1 = X1[1] * beta;
    }
    else {
      eta0 = csr_matrix_times_vector2(N[1], K, w_X0, v_X0, u_X0, beta);
      eta1 = csr_matrix_times_vector2(N[2], K, w_X1, v_X1, u_X1, beta);
    }
  }
  else {
    eta0 = rep_vector(0.0, N[1]);
    eta1 = rep_vector(0.0, N[2]);
  }
  if (has_intercept == 0 && dense_X) {
    real tmp = dot_product(xbar, beta);
    eta0 += tmp;
    eta1 += tmp;
  }
  if (has_offset == 1) {
    eta0 += offset0;
    eta1 += offset1;
  }
  if (K_smooth) {
    eta0 += S0 * beta_smooth;
    eta1 += S1 * beta_smooth;
  }
  if (special_case) for (i in 1:t) {
    eta0 += b[V0[i]];
    eta1 += b[V1[i]];
  }
  else if (t > 0) {
    eta0 += csr_matrix_times_vector2(N[1], q, w0, v0, u0, b);
    eta1 += csr_matrix_times_vector2(N[2], q, w1, v1, u1, b);
  }
