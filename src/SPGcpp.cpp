#include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec getEigenValues(arma::mat M) {
  return arma::eig_sym(M);
}

// [[Rcpp::export]]
arma::mat shrink_groupcpp(arma::mat C, arma::mat g_idx) {
  
  int V=g_idx.n_rows;
  int v;
  int n_rows = C.n_rows;
  int n_cols = C.n_cols;
  arma::mat R = mat(n_rows,n_cols);	
  arma::vec gnorm;
  arma::uvec g;
  
  for(v = 0; v < V; ++v){
    gnorm=sqrt(sum(pow(C.rows(g_idx(v,0)-1,g_idx(v,1)-1),2),0)).t();
    gnorm = clamp( gnorm, 1, R_PosInf);
    R.rows(g_idx(v,0)-1,g_idx(v,1)-1)=C.rows(g_idx(v,0)-1,g_idx(v,1)-1)/repmat(gnorm, g_idx(v,2), 1);
  }
  return R;
}

// [[Rcpp::export]]
double cal2normcpp(arma::mat A, arma::mat g_idx){
  
  int V = g_idx.n_rows;
  double s = 0;
  arma::mat B;
  arma::vec gnorm;
  for(int v = 0; v < V; ++v) {
    B = A.rows(g_idx(v,0)-1,g_idx(v,1)-1);
    gnorm=sqrt(sum(pow(B,2),0)).t();
    s = s+sum(gnorm);
  }
  return s;
}


// [[Rcpp::export]]
arma::vec soft_thresholding(arma::vec v, double lambdaL) {
  int len = v.size();
  arma::vec beta(len);
  for(int i = 0; i < len; ++i) {
    if(v(i) > lambdaL){
      beta(i) = v(i) - lambdaL;
    }
    else if(v(i) < -lambdaL){
      beta(i) = v(i) + lambdaL;
    }
    else{
      beta(i) = 0;
    }
  }
  return beta;
}


// [[Rcpp::export]]
List SPGcpp(arma::mat X, 
            arma::mat Y, 
            int J, 
            double gamma, 
            double lambda, 
            double Cnorm, 
            double mu,
            arma::sp_mat C,
            arma::vec w,
            arma::mat g_idx,
            double theta, 
            arma::vec beta,
            int maxiter,
            int display_iter,
            int N, 
            double tol){
  
  
  arma::mat XX = X.t()*X;
  arma::mat XY = X.t()*Y;
  int len = w.size();
  arma:: vec beta_new(len);
  arma::vec obj;
  arma::vec grad(len);
  double theta_new;
  double L;
  double CNorm = Cnorm;
  int iter;
  arma::mat Cmat = mat(C);
  
  // arma vec:: density 
  
  if (J < 10000) {
    L = max(getEigenValues(XX)) + pow(gamma,2)* CNorm / mu;
  }else{
    L = accu(pow(XX,2)) + pow(gamma,2) * CNorm / mu;
  }
  
  for(iter = 0; iter < maxiter; ++iter){
    
    // density(iter) = 1
    // Do the density portion later 
    // Add the line search part of the code 
    // Add the new grads and the new objective functions
    
    arma::mat A = shrink_groupcpp(Cmat * w / mu, g_idx);
    
    if (J < 2 * N && J < 10000){
      grad = XX * w - XY + Cmat.t()*A;
    }else{
      grad = X.t()*(X * w - Y) + Cmat.t()*A;
    }
    
    beta_new=soft_thresholding(w-(1/L)*grad, lambda/L);
    
    theta_new=(sqrt(pow(theta,4)+4*pow(theta,2))-pow(theta,2))/2;
    
    w=beta_new+((1-theta)/theta)*theta_new*(beta_new-beta);
    

    obj(iter)=sum(pow((Y-X*beta_new),2))/2+cal2normcpp(Cmat*beta_new, g_idx)+lambda*sum(abs(beta_new));
    
    if (iter > 10){
      if(abs(obj(iter)-obj(iter-1))/abs(obj(iter-1)) < tol) {
        break;
      }
    }

    if (iter==1 || iter % display_iter==0){
        Rcout<< "Iter, Objective";
    }
    beta = beta_new;
    theta = theta_new;
  }
  
  Rcout << "Iter, Objective";
  obj = obj.subvec(1,iter);
  
  return List::create( 
    _["obj"]  = obj, 
    _["iter"]  = iter, 
    _["beta"] = beta_new
  );
  
}