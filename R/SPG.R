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

  XX = Rfast::mat.mult(Rfast::transpose(X),X)
  XY = Rfast::mat.mult(X,Y)

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


# ## Algorithm in R LOOP
# 
# if(FALSE){
# 
# 
# 
#   XX = Matrix::t(X)%*%X
#   XY = Matrix::t(X)%*%Y
# 
#   # insert new calculation for lipschitz continuity
# 
#   if (J < 10000) {
#     L = eigs(A, 1, opts = list(retvec = FALSE))$values + gamma^2 * CNorm / mu
#   }else{
#     L = sum(XX^2) + gamma ^ 2 * CNorm / mu
#     rm(XX,XY)
#   }
# 
#   # tic
#   for(iter in 1:maxiter){
# 
#     ## density later
#     density[iter] = 1
# 
#     A = shrink_group(C %*% w / mu, g_idx)
# 
#     # include your new gradient functon here
# 
#     if (J < 2 * N && J < 10000){
#       grad = XX %*% w - XY + Matrix::t(C)%*%A;
#     }else{
#       grad = t(X)%*%(X %*% w - Y) + t(C)%*%A
#     }
# 
#     ## Final step to do
#     ## First do in R
#     ## Then do in Rcpp
#     beta_new=soft_threshold(w-(1/L)*grad, lambda/L);
#     beta_new=soft_thresholding(w-(1/L)*grad, lambda/L);
#     # density[iter]=density[iter]/J;
# 
#     theta_new=(sqrt(theta^4+4*theta^2)-theta^2)/2;
# 
#     w=beta_new+((1-theta)/theta)*theta_new*(beta_new-beta);
# 
#     ## include your new objective function here
# 
#     obj[iter]=sum((Y-X%*%beta_new)^2)/2+cal2norm(C%*%beta_new, g_idx)+lambda*sum(abs(beta_new));
# 
# 
# if (iter>10 & (abs(obj[iter]-obj[iter-1])/abs(obj[iter-1])<tol)) break;
# 
# 
#     if (verbose & (iter==1 | mod(iter, display_iter)==0)){
#       cat('Iter:',iter, ' Obj:',obj[iter], ' density:',density[iter],'\n')
#     }
# 
#     beta = beta_new
#     theta = theta_new
#   }
# 
#   # time = toc
#   if (verbose){
#     cat("\n", 'Iter:',iter, ' Obj:',obj[iter], "\n")
#   }
# 
#   obj = obj[1:iter]
#   density = density[1:iter]
#   return(list(obj = obj, density = density, beta = beta_new, iter = iter))
# }

