#include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double getmaxEigenValue(arma::mat M) {
  return max(arma::eig_sym(M));
}

// [[Rcpp::export]]
arma::mat shrink_groupcpp(arma::mat A, arma::mat g_idx) {
  
  int V=g_idx.n_rows;
  int v;
  int n_rows = A.n_rows;
  int n_cols = A.n_cols;
  arma::mat R = mat(n_rows,n_cols);	
  arma::vec gnorm;
  
  for(v = 0; v < V; ++v){
    gnorm=sqrt(sum(pow(A.rows(g_idx(v,0)-1,g_idx(v,1)-1),2),0)).t();
    if(gnorm(0,0) < 1){
      gnorm(0,0) = 1;
    }
    R.rows(g_idx(v,0)-1,g_idx(v,1)-1)=A.rows(g_idx(v,0)-1,g_idx(v,1)-1)/repmat(gnorm, g_idx(v,2), 1);
  }
  return R;
}

// [[Rcpp::export]]
double cal2normcpp(arma::mat A, arma::mat g_idx){
  
  int V = g_idx.n_rows;
  double s = 0;
  arma::vec gnorm;
  
  for(int v = 0; v < V; ++v) {
    gnorm=sqrt(sum(pow(A.rows(g_idx(v,0)-1,g_idx(v,1)-1),2),0)).t();
    s = s+sum(gnorm);
  }
  return(s);
}


// [[Rcpp::export]]
arma::mat soft_thresholdingcpp(arma::mat v, double lambdaL) {
  int len = v.n_rows;
  arma::mat beta(len,1);
  
  for(int i = 0; i < len; ++i) {
    if(v(i,0) > lambdaL){
      beta(i,0) = v(i,0) - lambdaL;
    }
    else if(v(i,0) < -lambdaL){
      beta(i,0) = v(i,0) + lambdaL;
    }
    else{
      beta(i,0) = 0;
    }
  }
  return beta;
}


// [[Rcpp::export]]
List SPGcpp(arma::mat X, 
            arma::mat Y, 
            arma:: mat XX,
            arma:: mat XY,
            int J, 
            double gamma, 
            double lambda, 
            double CNorm, 
            double mu,
            arma::sp_mat C,
            arma::mat w,
            arma::mat g_idx,
            double L,
            double theta, 
            arma::mat beta,
            int maxiter,
            int display_iter,
            int N, 
            double tol){
  
  
  int len = w.n_rows;
  arma::mat beta_new(len,1);
  arma::mat obj(maxiter,1);
  arma::mat grad(len,1);
  double theta_new;
  int iter;
  arma::mat Cmat = mat(C);
  
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

    beta_new = soft_thresholdingcpp(w-grad/L, lambda/L);
      

    theta_new = (sqrt(pow(theta,4)+4*pow(theta,2))-pow(theta,2))/2;

    w = beta_new + ((1-theta)/theta)*theta_new*(beta_new-beta);

    obj.row(iter) =  sum(pow(Y - X*beta_new,2))/2 + cal2normcpp(Cmat*beta_new, g_idx)+lambda*sum(abs(beta_new));

    if (iter > 10){
      if( abs(obj(iter,0)-obj(iter-1,0))/abs(obj(iter-1,0)) < tol  ) {
        // break;
      }
    }

    if (iter % display_iter == 0){
      Rprintf("Iter: %i, Objective %f \n", iter, obj(iter,0));
      }

    beta = beta_new;
    theta = theta_new;
  }
  iter = iter-1;
  Rprintf("Iter: %i, Objective %f \n", iter, obj(iter,0));
  obj = obj.col(0); //.subvec(0,iter);

  return List::create(
    _["obj"]  = obj,
    _["iter"]  = iter,
    _["beta"] = beta_new
  );
  
}