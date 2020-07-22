# this code assumes that variables in cov.mat are ordered as (pa.x,var.x,others)
derivative.of.the.modified.cov.mat <- function(cov.mat,q1){
   
  nvar <- nrow(cov.mat)
  if (q1 == 0) return(diag(nvar*(nvar+1)/2))
  
  pa <- 1:q1
  z <- 1:(q1+1)
  w <- (q1+2):nvar
  
  cov.mat.zz.inv <- solve(cov.mat[z,z])
  cov.mat.pa.inv <- solve(cov.mat[pa,pa])
  modified.cov.mat.zz <- cov.mat[z,z]
  modified.cov.mat.zz[q1+1,pa] <- rep(0,q1)
  modified.cov.mat.zz[pa,q1+1] <- rep(0,q1)
  modified.cov.mat.zz[q1+1,q1+1] <- cov.mat[q1+1,q1+1] - cov.mat[q1+1,pa]%*%cov.mat.pa.inv%*%cov.mat[pa,q1+1]
  
  H <- modified.cov.mat.zz%*%cov.mat.zz.inv
  R <- cov.mat[w,z]%*%cov.mat.zz.inv
  Q <- rbind(H,R%*%(H - diag(q1+1)))
  Q <- cbind(Q,rbind(matrix(0,nrow=q1+1,ncol=nvar-q1-1),diag(nvar-q1-1)))
  
  derivative.of.modified.cov.mat.at.idx <- function(idx){
    r <- idx[1]
    s <- idx[2]
    E <- matrix(0,nrow=nvar,ncol=nvar)
    E[r,s] <- 1
    E[s,r] <- 1
    
    if (min(r,s)>q1+1){
      return(sm2vec(Q%*%E,diag=T))
    }else{#compute derivative of Q
      derivative.of.Q <- matrix(0,nrow=nvar,ncol=nvar)
      if (max(r,s)>q1+1){# derivative of H is 0
        derivative.of.Q[w,z] <- E[w,z]%*%cov.mat.zz.inv%*%(H - diag(q1+1))
      }else{
        #first compute derivative of H
        derivative.of.H  <- matrix(0,nrow=q1+1,ncol=q1+1)
        #to this end, compute derivative.of.modified.cov.mat.zz
        derivative.of.modified.cov.mat.zz <- matrix(0,nrow=q1+1,ncol=q1+1)
        if (max(r,s)<=q1){
          derivative.of.modified.cov.mat.zz[pa,pa] <- E[pa,pa]
          derivative.of.modified.cov.mat.zz[q1+1,q1+1] <- cov.mat[q1+1,pa]%*%cov.mat.pa.inv%*%E[pa,pa]%*%cov.mat.pa.inv%*%cov.mat[pa,q1+1]
        }else{
          derivative.of.modified.cov.mat.zz[q1+1,q1+1] <- if(r==q1+1 & s==q1+1) 1 else (-2*cov.mat.pa.inv[min(r,s),]%*%cov.mat[pa,q1+1])
        }# computation of derivative.of.modified.cov.mat.zz ends here
        derivative.of.H <- derivative.of.modified.cov.mat.zz%*%cov.mat.zz.inv - modified.cov.mat.zz%*%cov.mat.zz.inv%*%E[z,z]%*%cov.mat.zz.inv
        derivative.of.Q[z,z] <- derivative.of.H
        derivative.of.Q[w,z] <- - R%*%E[z,z]%*%cov.mat.zz.inv%*%( H - diag(q1+1)) + R%*%derivative.of.H
      }# computation of deviative of Q ends here
      return(sm2vec(derivative.of.Q%*%cov.mat + Q%*%E,diag=T))
    }
  }

  halfvec.idx <- sm.index(cov.mat,diag=TRUE)
  halfvec.idx <- halfvec.idx[,c(2,1)]
  return(apply(halfvec.idx,1,derivative.of.modified.cov.mat.at.idx))
}