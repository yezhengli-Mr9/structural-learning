do.intervention <- function(cov.mat,pa2,ind2){
  if (length(pa2)>0){  # must cut edges into x2
    est.B <- t(solve(chol(cov.mat)))
    est.B[ind2,pa2] <- 0
    #cov.mat2 <- solve(t(est.B2)%*%est.B2)
    cov.mat <- solve(est.B)%*%t(solve(est.B))
  }
  return(cov.mat)
  
}
calculate.cov.mat.of.modified.halfvec.cov.mat <- function(cov.mat,ivar,ipaset,npa,cov.of.halfvec.cov.mat=NULL){

  nvar <- NROW(cov.mat)
  all.var <- 1:nvar
  ind <- 1:nvar
  ivar.k <- ivar
  ipaset.k <- ipaset
  #indy.k <- indy
  cov.mat.k<-cov.mat
  if (is.null(cov.of.halfvec.cov.mat))
  cov.of.halfvec.cov.mat <- calculate.cov.of.cov.mat(cov.mat)
  #rownames(cov.of.halfvec.cov.mat)<-1:NROW(cov.of.halfvec.cov.mat)
  cov.of.halfvec.modified.cov.mat <- cov.of.halfvec.cov.mat
  halfvec.idx <- sm.index(cov.mat,diag=TRUE)
  halfvec.idx<-halfvec.idx[,c(2,1)]
  vec.idx <- cbind(rep(1:nvar,nvar),rep(1:nvar,rep(nvar,nvar)))
  transform.vec2vech <- Matrix(0,nrow=NROW(halfvec.idx),ncol=NROW(vec.idx),sparse=TRUE)
  for (i in 1:NROW(halfvec.idx)){
    transform.vec2vech[i,which(vec.idx[,1]%in%halfvec.idx[i,1] & vec.idx[,2]%in%halfvec.idx[i,2])] <- 1
    #col.idx <- c(col.idx,which(vec.idx[,1]%in%halfvec.idx[i,1] & vec.idx[,2]%in%halfvec.idx[i,2]))
  }
  transform.vech2vec <- Matrix(0,nrow=NROW(vec.idx),ncol=NROW(halfvec.idx),sparse=TRUE)
  for (i in 1:NROW(vec.idx)){
    transform.vech2vec[i,which((halfvec.idx[,1]%in%vec.idx[i,1] & halfvec.idx[,2]%in%vec.idx[i,2]) |(halfvec.idx[,1]%in%vec.idx[i,2] & halfvec.idx[,2]%in%vec.idx[i,1]) )] <- 1
    #col.idx <- c(col.idx,which(vec.idx[,1]%in%halfvec.idx[i,1] & vec.idx[,2]%in%halfvec.idx[i,2]))
  }
  
  for (k in 1:length(ivar)){
    if (npa[k]>0){
      indk <- ivar.k[k]
      if (k==1){
        pak <- ipaset.k[1:npa[1]]
      }else{
          pak <- ipaset.k[(sum(npa[1:(k-1)])+1):sum(npa[1:k])]
      }
      
      ind <- c(pak,indk,setdiff(ind,c(pak,indk)))
      pak <- match(pak,ind)
      indk <- match(indk,ind)
      cov.mat.k <- cov.mat.k[ind,ind]
      all.var <- all.var[ind]
      ipaset.k <- match(ipaset.k,ind)
      ivar.k <- match(ivar.k,ind)
      
      perm.mat <- as(ind,"pMatrix")
      halfvec.perm.mat <- transform.vec2vech%*%kronecker(perm.mat,as.matrix(perm.mat))%*%transform.vech2vec
      cov.of.halfvec.modified.cov.mat <- halfvec.perm.mat%*%cov.of.halfvec.modified.cov.mat%*%t(halfvec.perm.mat)
      gmat<- derivative.of.the.modified.cov.mat(cov.mat.k,indk,pak,halfvec.idx)
      # gmat<- transform.mat%*%gmat%*%t(transform.mat)
      cov.of.halfvec.modified.cov.mat <- t(gmat)%*%cov.of.halfvec.modified.cov.mat%*%gmat
      cov.mat.k <- do.intervention(cov.mat.k,pak,indk)
    }
  }
  perm.mat <- as(order(all.var),"pMatrix")
  modified.cov.mat <- perm.mat%*%cov.mat.k%*%t(perm.mat)
  halfvec.perm.mat <- transform.vec2vech%*%kronecker(perm.mat,perm.mat)%*%transform.vech2vec
  cov.of.halfvec.modified.cov.mat <- halfvec.perm.mat%*%cov.of.halfvec.modified.cov.mat%*%t(halfvec.perm.mat)
return(list(cov.of.halfvec.modified.cov.mat=cov.of.halfvec.modified.cov.mat,modified.cov.mat=modified.cov.mat))
 
}