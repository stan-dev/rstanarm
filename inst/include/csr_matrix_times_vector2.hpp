#ifndef RSTANARM__CSR_MATRIX_TIMES_VECTOR2_HPP
#define RSTANARM__CSR_MATRIX_TIMES_VECTOR2_HPP

/*
 * This works exactly like csr_matrix_times_vector but faster and less safe
 */
template <typename T1, typename T2>
inline
Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
              Eigen::Dynamic, 1>
csr_matrix_times_vector2(const int& m,
                         const int& n,
                         const Eigen::Matrix<T1, Eigen::Dynamic, 1>& w,
                         const std::vector<int>& v,
                         const std::vector<int>& u,
                         const Eigen::Matrix<T2, Eigen::Dynamic, 1>& b,
                         std::ostream* pstream__) {
  Eigen::Map<const Eigen::SparseMatrix<T1,Eigen::RowMajor> >
    sm(m, n, w.size(), &u[0], &v[0], &w[0]);
  return sm * b;
}

/* This specialization is slower than the above templated version
stan::math::vector_v
csr_matrix_times_vector2(const int& m,
                         const int& n,
                         const Eigen::VectorXd& w,
                         const std::vector<int>& v,
                         const std::vector<int>& u,
                         const stan::math::vector_v& b,
                         std::ostream* pstream__) {
  Eigen::Map<const Eigen::SparseMatrix<double,Eigen::RowMajor> >
    sm(m, n, w.size(), &u[0], &v[0], &w[0]);
  Eigen::VectorXd result = sm * value_of(b);
  stan::math::vector_v out(m);
  int pos = 0;
  for (int k = 0; k < m; k++) {
    int nnz = u[k + 1] - u[k];
    if (nnz > 0) {
      std::vector<stan::math::var> operands(nnz);
      std::vector<double> gradients(nnz);
      for (int j = 0; j < nnz; j++) {
        operands[j] = b.coeff(v[pos]);
        gradients[j] = w.coeff(pos++);
      }
      out.coeffRef(k) = precomputed_gradients(result.coeff(k),
                                              operands, gradients);
    }
    else out.coeffRef(k) = 0.0;
  }
  return out;
}
*/
#endif
