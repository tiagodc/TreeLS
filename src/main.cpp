// #include <Rcpp.h>
// // #include <RcppEigen.h>
//
// // // [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::plugins("cpp11")]]
//
// using namespace Rcpp;
// using Eigen::Map;                       // 'maps' rather than copies
// using Eigen::MatrixXd;                  // variable size matrix, double precision
// using Eigen::VectorXd;                  // variable size vector, double precision
// using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
//
// NumericVector timesTwo(NumericVector x) {
//   return x * 2;
// }
//
// VectorXd getEigenValues(Map<MatrixXd> M) {
//   SelfAdjointEigenSolver<MatrixXd> es(M);
//   return es.eigenvalues();
// }
