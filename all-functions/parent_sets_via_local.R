parent.sets.via.local <- function(x,amat){ #x must be a scaler
 # amat[which(amat != 0)] <- 1
  amat.undir <- amat*t(amat) 
  amat.dir <- amat - amat.undir 
  pa.dir <- which(amat.dir[,x]!=0)
  pa <- which(amat.undir[,x] != 0)
  if (length(pa)==0)
    return(list(pa.dir))
  paset <- list(pa.dir)
  if (length(pa)==1){
    paset <- c(paset,list(c(pa,pa.dir))) 
  }else{
    for (i in 1:length(pa)){
      pa.tmp <- combn(pa, i, simplify = TRUE)
      n.comb <- ncol(pa.tmp)
      for (j in 1:n.comb) {
        pa.t <- pa.tmp[, j]
        new.coll <- FALSE
        if (length(pa.t)>1){
          tmp <- amat.undir[pa.t,pa.t]
          diag(tmp) <- 1
          new.coll <- (min(tmp)==0)
        }
        if (!new.coll) paset <- c(paset,list(c(pa.t,pa.dir)))
      }
    }
  }
  
  return(paset)
}