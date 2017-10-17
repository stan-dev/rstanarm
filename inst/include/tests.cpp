// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <stan/math.hpp> 
#include "meta_header.hpp"

// [[Rcpp::export]]
Eigen::VectorXd
csr_matrix_times_vector2_test(const int& m,
                              const int& n,
                              const Eigen::VectorXd& w,
                              const std::vector<int>& v,
                              const std::vector<int>& u,
                              const Eigen::VectorXd& b) {
  return csr_matrix_times_vector2(m,n,w,v,u,b,0);
}
