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
