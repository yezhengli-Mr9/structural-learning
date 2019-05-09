asymp.var.MCD <- function(cov.mat,intervention.set,var.y,pasets,cov.of.halfvec.cov.mat=NULL){
  
  
  p <- nrow(cov.mat)
  imp.var <- unique(c(intervention.set,unlist(pasets),var.y))
  if (length(imp.var) < p){
    cov.mat <- cov.mat[imp.var,imp.var]
    intervention.set <- match(intervention.set,imp.var)
    var.y <- match(var.y,imp.var)
    pasets <- lapply(pasets,function(x) match(unlist(x),imp.var))
  }
  
  
  if (is.null(cov.of.halfvec.cov.mat) | length(imp.var) < p)
    cov.of.halfvec.cov.mat <- calculate.cov.of.cov.mat(cov.mat)
  
  p <-  length(imp.var)
  halfvec.idx <- sm.index(cov.mat,diag=TRUE)
  halfvec.idx <- halfvec.idx[,c(2,1)]
  vec.idx <- as.matrix(expand.grid(1:p,1:p))
  transform.vec2vech <- Matrix(0,nrow=NROW(halfvec.idx),ncol=NROW(vec.idx),sparse=TRUE)
  transform.vech2vec <- Matrix(0,nrow=NROW(vec.idx),ncol=NROW(halfvec.idx),sparse=TRUE)
  for (i in 1:NROW(halfvec.idx))
    transform.vec2vech[i,which(vec.idx[,1]==halfvec.idx[i,1] & vec.idx[,2]==halfvec.idx[i,2])] <- 1
  for (i in 1:NROW(vec.idx))
    transform.vech2vec[i,which(halfvec.idx[,1]==max(vec.idx[i,]) & halfvec.idx[,2]==min(vec.idx[i,]))] <- 1
  
  asymp.var.mat.perm.fun <- function(ind,asymp.var.mat){
    perm.mat <- as(ind,"pMatrix")
    halfvec.perm.mat <- transform.vec2vech%*%kronecker(perm.mat,perm.mat)%*%transform.vech2vec
    return(halfvec.perm.mat%*%asymp.var.mat%*%t(halfvec.perm.mat))
  }
  
  compute.asymp.var.mat <- function(cov.mat.temp,asymp.var.mat,var.x.temp){#compute the covariance matrix of the transformed linear SEM
    #Obtain the parents set of var.x assuming that if var.x does not belong to intervention.set then var.x has no parents
    pa.x <-  if (is.element(var.x.temp,intervention.set)) unlist(pasets[match(var.x.temp,intervention.set)])
    if (length(pa.x)==0) return(asymp.var.mat)
    # order the variables in cov.mat as necessary
    ind <- c(pa.x,var.x.temp,(1:nrow(cov.mat.temp))[-c(pa.x,var.x.temp)]) 
    cov.mat.temp <- cov.mat.temp[ind,ind]
    Lambda <- derivative.of.the.modified.cov.mat(cov.mat.temp,length(pa.x))
    asymp.var.mat <- Lambda%*%asymp.var.mat.perm.fun(ind,asymp.var.mat)%*%t(Lambda)
    return(asymp.var.mat.perm.fun(order(ind),asymp.var.mat))
  }
  
  Gamma.mat <- cov.of.halfvec.cov.mat
  modified.cov.mat <- cov.mat
  
  for (i in 1:length(intervention.set)){
    Gamma.mat <- compute.asymp.var.mat(modified.cov.mat,Gamma.mat,intervention.set[i])
    modified.cov.mat <- MCD(modified.cov.mat,intervention.set[i],var.y,pasets[i],return.modified.cov.mat=TRUE)
  }
  
  H <- Matrix(0,nrow=length(intervention.set),ncol=p*(p+1)/2,sparse=T)
  for (i in 1:length(intervention.set)){
    H[i,which(halfvec.idx[,1]==max(intervention.set[i],var.y) & halfvec.idx[,2]==min(intervention.set[i],var.y))] <- 1/modified.cov.mat[intervention.set[i],intervention.set[i]]
    H[i,which(halfvec.idx[,1]==intervention.set[i] & halfvec.idx[,2]==intervention.set[i])] <- -(modified.cov.mat[intervention.set[i],var.y])/(modified.cov.mat[intervention.set[i],intervention.set[i]])^2
  }
  
  asymp.var.mat <- as.matrix(H%*%Gamma.mat%*%t(H))
  
  for (i in 1:length(intervention.set)){
    if (is.element(var.y,unlist(pasets[i]))){
      asymp.var.mat[,intervention.set[i]] <- 0; asymp.var.mat[intervention.set[i],] <- 0;
    }
  }
  return(asymp.var.mat)  
}
