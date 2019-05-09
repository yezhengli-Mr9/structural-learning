asymp.var.RR <- function(cov.mat,intervention.set,var.y,pasets,cov.of.halfvec.cov.mat=NULL){
  
  if (length(intervention.set)>2)
    cat("\n This function is only for at max 2 interventions \n")
  
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
  
  diff.fun <- function(x,y){
    derivative.mat <- matrix(0,nrow=p,ncol=p)
    z <- c(unlist(pasets[match(x,intervention.set)]),x)
    nz <- length(z)
    if (!is.element(y,z)){
      cov.mat.zz.inv <- solve(cov.mat[z,z])
      H <- cov.mat.zz.inv%*%cov.mat[z,y]%*%cov.mat.zz.inv[nz,]
      derivative.mat[z,z] <- -(H + t(H) - diag(diag(matrix(H,nrow=nz)),nrow=nz))
      derivative.mat[z,y] <- cov.mat.zz.inv[,nz]
      derivative.mat[y,z] <- cov.mat.zz.inv[nz,]
    }
    return(derivative.mat)
  }
  
  
  if (length(intervention.set) == 1){
    Lambda  <- sm2vec(diff.fun(intervention.set[1],var.y),diag=TRUE)
    return(matrix(Lambda,nrow=1)%*%cov.of.halfvec.cov.mat%*%matrix(Lambda,ncol=1))
  }else{
     Lambda <- rbind(sm2vec(diff.fun(intervention.set[1],var.y),diag=TRUE),sm2vec(diff.fun(intervention.set[2],var.y),diag=TRUE),
                    sm2vec(diff.fun(intervention.set[1],intervention.set[2]),diag=TRUE),sm2vec(diff.fun(intervention.set[2],intervention.set[1]),diag=TRUE))

     #define function for calculating single intervention effect of var.x on var.y
     RR.single.intervention <- function(var.x,var.y){
       if (var.x == var.y) return(0)
       #Obtain the parents set of var.x assuming that if var.x does not belong to intervention.set then var.x has no parents
       pa.x <-  if (is.element(var.x,intervention.set)) unlist(pasets[match(var.x,intervention.set)]) 
       if(is.element(var.y,pa.x)) return(0)
       #the single intervention effect is the regression coefficient of var.x in the regression var.y ~ var.x + pa(var.x)
       return(solve(cov.mat[c(var.x,pa.x),c(var.x,pa.x)])[1,]%*%cov.mat[c(var.x,pa.x),var.y])
     }
     
    causal.effects <- c(RR.single.intervention(intervention.set[1],var.y),RR.single.intervention(intervention.set[2],var.y),
                            RR.single.intervention(intervention.set[1],intervention.set[2]),RR.single.intervention(intervention.set[2],intervention.set[1]))
    
  Fmat <- matrix(c(1,-causal.effects[4],-causal.effects[3],1,-causal.effects[2],0,0,-causal.effects[1]),nrow=2)
     asymp.var <- Fmat%*%Lambda%*%cov.of.halfvec.cov.mat%*%t(Lambda)%*%t(Fmat)
     return(asymp.var)
  }
                              
}

