  vector[N[1]] eta0;
  vector[N[2]] eta1;
  if (K > 0) {
    eta0 <- X0 * beta;
    eta1 <- X1 * beta;
  }
  else {
    eta0 <- rep_vector(0.0, N[1]);
    eta1 <- rep_vector(0.0, N[2]);
  }
  if (has_intercept == 0) {
    real tmp;
    tmp <- dot_product(xbar, beta);
    eta0 <- eta0 + tmp;
    eta1 <- eta1 + tmp;
  }
  if (has_offset == 1) {
    eta0 <- eta0 + offset0;
    eta1 <- eta1 + offset1;
  }
  if (t > 0) {
    eta0 <- eta0 + csr_matrix_times_vector(N[1], q, w0, v0, u0, b);
    eta1 <- eta1 + csr_matrix_times_vector(N[2], q, w1, v1, u1, b);
  }
