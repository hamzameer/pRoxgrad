gentoy_group <- function(n, ng, g_size, overlap){

  #### Return the following elements
  # X: design matrix
  # Y: output
  # Tm: ng \times p, indicate each group contains which variables
  # Tw: weight for each group
  # w: true regression coefficients

  ### With the following inputs
  # ng  # number of groups
  # g_size  # group size
  # overlap  # number of variables overlapped between two consecutive groups
  # n     # sample size
  # p=ng*(g_size-overlap)+overlap   # total number of variables # 4510

  d=ng*(g_size-overlap)+overlap

  r_T=kronecker(1:ng, rep(1,g_size))
  c_T=rep(0,ng*g_size)
  s=1
  for (g in 1:ng){
    c_T[((g-1)*g_size+1):(g*g_size)]=s:(s+g_size-1)
    s=s+g_size-overlap
  }

  Tm = Matrix::sparseMatrix(i=r_T, j=c_T, x= 1, dims=c(ng, d)) # ng x d | 50 x 4510

  Tw=rep(1,ng)  # uniform weight

  tmp=1:d
  w= sapply(tmp, function(x) ((-1)^x)*(exp(-(x-1)/100)) )
  sn2 = 1;  # signal to noise ratio

  X = matrix(rnorm(n*d), nrow = n)
  Y=  X%*%w + t(sn2%*%rnorm(n))
  return(list(X=X,Y=Y,Tm=Tm,Tw=Tw,w=w))

}
