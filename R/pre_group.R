pre_group <- function(Tm,Tw){

  # Return the values of
  # C | Matrix C
  # g_idx |
  # TauNorm | Square of ||C||

  V = dim(Tm)[1]; K = dim(Tm)[2]
  sum_col_T=Matrix::rowSums(Tm)
  SV=sum(sum_col_T)
  csum=cumsum(sum_col_T)

  g_idx = cbind(seq(1,csum[49]+1,100), csum, sum_col_T)
  # each row is the range of the group

  J = rep(0,SV)
  W =  rep(0,SV)

  for(v in 1:V){
    J[g_idx[v,1]:g_idx[v,2]]=pracma::finds(Tm[v,])
    W[g_idx[v,1]:g_idx[v,2]]=Tw[v]
  }

  C=Matrix::sparseMatrix(i = 1:SV, j = J, x = W, dims = c(SV, K))

  TauNorm=Matrix::sparseMatrix(i=1:V,j=1:V, x = Tw, dims = c(V,V))%*%Tm
  TauNorm= max(Matrix::colSums(TauNorm^2))
  return(list(C = C, g_idx = g_idx, Cnorm = TauNorm))
}
