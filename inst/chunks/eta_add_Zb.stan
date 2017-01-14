    if (special_case) for (i in 1:t) eta = eta + b[V[i]];
    else eta = eta + csr_matrix_times_vector(N, q, w, v, u, b);
