#include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double cal2normcpp(arma::sp_mat A, NumericMatrix g_idx){
  
  int V = A.n_rows;
  int C = A.n_cols;
  double s = 0;
  IntegerVector idx;
  NumericMatrix B;
  for(int i = 0; i < V; ++i) {
    idx = seq(g_idx(i,1),g_idx(i,2));
    // B = A(Range(g_idx(i,1)-1, C-1), Range(g_idx(i,2)-1, C-1));
    B = A.row(idx);
    int n = idx.size()*C;
    double gnorm = 0;
    for(int j = 0; j < n; ++j){
      gnorm += pow(B[j],2.0);
    }
    gnorm=std::sqrt(gnorm);
    s+=gnorm;
  }
  return(s);
}