#ifndef RSTANARM__CODOLS_HPP
#define RSTANARM__CODOLS_HPP

/*
 * Compute ordinary least squares coefficients,
 * even in the situation where X is rank deficient
 * See https://eigen.tuxfamily.org/dox/classEigen_1_1CompleteOrthogonalDecomposition.html
 */

template <typename T2__, typename T3__>
inline
Eigen::Matrix<typename boost::math::tools::promote_args<T2__, T3__>::type,
              Eigen::Dynamic, 1>
CODOLS(const Eigen::Matrix<T2__, Eigen::Dynamic, Eigen::Dynamic>& X,
       const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& y,
       std::ostream* pstream__) {
  typename boost::math::tools::promote_args<T2__, T3__>::type T1__;
  using namespace Eigen;
  CompleteOrthogonalDecomposition<MatrixXd> cod(X);
  return cod.solve(y);
}

inline
auto
CODOLS(const Eigen::Matrix<double, -1, -1, 0, -1, -1>& X,
       const Eigen::Map<Eigen::Matrix<double, -1, 1, 0, -1, 1>, 0, Eigen::Stride<0, 0>>& y,
       std::ostream* pstream__) {
  using namespace Eigen;
  CompleteOrthogonalDecomposition<MatrixXd> cod(X);
  return cod.solve(y).eval();
}
#endif
