    if (special_case) eta = eta + b[v];
    else eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);
