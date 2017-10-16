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
  typedef typename boost::math::tools::promote_args<T1, T2>::type result_t;
  Eigen::Map<const Eigen::SparseMatrix<T1,Eigen::RowMajor> >
    sm(m, n, w.size(), &u[0], &v[0], &w[0]);
  return sm * b;
}
