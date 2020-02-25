soft_threshold <- function(v, lambda_L){
  sapply(v,function(x) sign(x)*max(0,abs(x) - lambda_L))
}