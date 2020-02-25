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
