calculate.cov.of.cov.mat <- function(cov.mat){
  p <- nrow(cov.mat)
  halfvec.idx <- sm.index(cov.mat,diag=TRUE)
  halfvec.idx <- halfvec.idx[,c(2,1)]
  vec.idx <- as.matrix(expand.grid(1:p,1:p))
  transform.vec2vech <- Matrix(0,nrow=NROW(halfvec.idx),ncol=NROW(vec.idx),sparse=TRUE)
  for (i in 1:NROW(halfvec.idx))
    transform.vec2vech[i,which(vec.idx[,1]==halfvec.idx[i,1] & vec.idx[,2]==halfvec.idx[i,2])] <- 1
  rowidx <- apply(vec.idx,1,function(x) p*(x[1]-1)+x[2]) 
  colidx <- apply(vec.idx,1,function(x) p*(x[2]-1)+x[1])
  comm.mat <- spMatrix(nrow=p^2,ncol=p^2,i=rowidx,j=colidx,x=rep(1,p^2))
  return(as.matrix(transform.vec2vech%*%(Diagonal(p^2)+comm.mat)%*%kronecker(cov.mat,cov.mat)%*%t(transform.vec2vech)))
}


