#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector soft_thresholding(NumericVector v, double lambdaL) {
  
  int len = v.size();
  NumericVector beta(len);
  for(int i = 0; i < len; ++i) {
    if(v[i] > lambdaL){
      beta[i] = v[i] - lambdaL;
    }
    else if(v[i] < -lambdaL){
      beta[i] = v[i] + lambdaL;
    }
    else{
      beta[i] = 0;
    }
  }
  return beta;
}