MCD <- function(cov.mat,intervention.set,var.y,pasets,return.modified.cov.mat=FALSE){
  if (is.element(var.y,intervention.set) & !return.modified.cov.mat) return(rep(0,length(intervention.set)))
  if (length(intervention.set)==1 & is.element(var.y,unlist(pasets)) & !return.modified.cov.mat) return(0)
  
  if (!return.modified.cov.mat){
    imp.var <- unique(c(intervention.set,unlist(pasets),var.y))
    if (length(imp.var) < nrow(cov.mat)){
      cov.mat <- cov.mat[imp.var,imp.var]
      intervention.set <- match(intervention.set,imp.var)
      var.y <- match(var.y,imp.var)
      pasets <- if(length(intervention.set)>1) lapply(pasets,function(x) match(unlist(x),imp.var)) else match(unlist(pasets),imp.var)
    }
  }

  do.Cholesky.modification <- function(x){
    pa.x <-  if(length(intervention.set)>1) pasets[[match(x,intervention.set)]] else unlist(pasets)
    if (length(pa.x)==0) return(cov.mat)
    ind <- c(pa.x,x,(1:nrow(cov.mat))[-c(pa.x,x)]) 
    cov.mat <- cov.mat[ind,ind]
    x <- match(x,ind)
    pa.x <- match(pa.x,ind)
    temp <- gchol(cov.mat)
    Lower.tri.mat <- solve(as.matrix(temp))
    Diag.mat <- diag(diag(temp))
    Lower.tri.mat[x,pa.x] <- 0
    cov.mat <- solve(Lower.tri.mat)%*%Diag.mat%*%t(solve(Lower.tri.mat))
    return(cov.mat[order(ind),order(ind)]) 
  }

    for (i in 1:length(intervention.set))
      cov.mat <- do.Cholesky.modification(intervention.set[i])
  
  if(return.modified.cov.mat) return(cov.mat)
  MCD.estimate <- function(x) if (is.element(var.y,unlist(pasets[match(x,intervention.set)]))) return(0) else return(cov.mat[var.y,x]/cov.mat[x,x])
  return(sapply(intervention.set,MCD.estimate))
  
}