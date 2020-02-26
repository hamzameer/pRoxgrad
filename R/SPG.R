SPG <- function(prob = 'group',Y, X, gamma, lambda, C, CNorm, option  = NULL, g_idx) {
  # Smoothing Proximal Gradient based on FISTA for medium-scale problem
  # with a pre-computed Lipschitz constant (without line-search)

  # problem to solve either 'group'
  # X inputs
  # gamma:  regularization for group norm
  # lambda: regularization for L1 penalty
  # C:  \sum_|g| by J matrix or |E| by J matrix
  # CNorm: ||C||^2
  # option: maxiter; tol, b0 ;
  # g_idx: n_group by 2, group index

  N = dim(X)[1]; J = dim(X)[2]

  if (!is.null(option$maxiter)){
    maxiter = option$maxiter
  }else maxiter = 10000

  if (!is.null(option$verbose)){
    verbose = option$verbose
  }else verbose = TRUE

  display_iter = 10

  if (!is.null(option$tol)){
    tol = option.tol
  }else tol = 1e-9

  if (!is.null(option$b_init)){
    beta = option$b_init
  }else beta = matrix(0, nrow = J, ncol=1)

  if(!is.matrix(beta)){
    beta = as.matrix(beta, ncol = 1)
  }

  if (!is.null(option$mu)){
    mu = option$mu
  }else mu = 1e-3

  C = C * gamma

  w = beta

  theta = 1

  XX = crossprod(X)
  XY = crossprod(X,Y)

  if (J < 10000) {
    L = RSpectra::eigs(XX, 1, opts = list(retvec = FALSE))$values + gamma^2 * CNorm / mu
  }else{
    L = sum(XX^2) + gamma ^ 2 * CNorm / mu
    rm(XX,XY)
  }
  
  solve <- pRoxgrad::SPGcpp(X, Y, XX, XY, J, gamma, lambda, CNorm, mu, C, w, g_idx, L,
                  theta, beta, maxiter=500, display_iter=10, N, tol)
  return (solve)
}

SPGR <- function(X, Y, XX, XY,
                 J, gamma, lambda,
                 CNorm, mu, C, w,
                 g_idx, L,
                 theta, beta,
                 maxiter=500,
                 display_iter=10,
                 N, tol) {
  
  obj = rep(0,maxiter)
  
  
  for(iter in 1:maxiter){
    
    
    A = shrink_group(C %*% w / mu, g_idx)
    
    if (J < 2 * N && J < 10000){
      grad = XX %*% w - XY + Matrix::t(C)%*%A;
    }else{
      grad = t(X)%*%(X %*% w - Y) + Matrix::t(C)%*%A
    }
    
    # beta_new=soft_threshold(w-(1/L)*grad, lambda/L);
    
    beta_new=pRoxgrad::soft_thresholdingcpp(as.matrix(w-(1/L)*grad), lambda/L);
    beta_new = as.vector(beta_new)
    

    theta_new=(sqrt(theta^4+4*theta^2)-theta^2)/2;
    
    w=beta_new+((1-theta)/theta)*theta_new*(beta_new-beta);
    

    obj[iter]=sum((Y-X%*%beta_new)^2)/2+cal2norm(C%*%beta_new, g_idx)+lambda*sum(abs(beta_new));
    
    # if (iter>10 & (abs(obj[iter]-obj[iter-1])/abs(obj[iter-1])<tol)) break;
    
    
    if (iter==1 | pracma::mod(iter, display_iter)==0){
      cat('Iter:',iter, ' Obj:',obj[iter],'\n')
    }
    
    
    beta = beta_new
    theta = theta_new
  }
  
  
  # time = toc
  cat("\n", 'Iter:',iter, ' Obj:',obj[iter], "\n")
  
  obj = obj[1:iter]

  return(list(obj = obj, beta = beta_new, iter = iter))
}
  

shrink_group <- function(A, g_idx){
  
  V=nrow(g_idx);
  R=matrix(0, dim(A)[1], dim(A)[2])
  
  for(v in 1:V){
    idx=g_idx[v,1]:g_idx[v,2]
    gnorm=sqrt(sum(A[idx,]^2))
    gnorm[gnorm<1]=1;
    R[idx,]=A[idx,]/pracma::repmat(gnorm,  g_idx[v,3], 1);
  }
  return(R)
}

soft_threshold <- function(v, lambda_L){
  sapply(v,function(x) sign(x)*max(0,abs(x) - lambda_L))
}

cal2norm <- function(A, g_idx) {
  V=nrow(g_idx)
  s=0
  for(v in 1:V){
    idx=g_idx[v,1]:g_idx[v,2];
    gnorm=sqrt(sum(A[idx,]^2));
    s=s+sum(gnorm);
  }
  return(s)
}






