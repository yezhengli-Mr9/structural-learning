RR <- function(cov.mat,intervention.set,var.y,pasets){
  
   adjusted.regression <- function(x,y){
      if (x == y) return(0)
      pa.x <-  if(length(intervention.set)>1) pasets[[match(x,intervention.set)]] else unlist(pasets)
      if(is.element(y,pa.x)) return(0)
      return(solve(cov.mat[c(x,pa.x),c(x,pa.x)],cov.mat[c(x,pa.x),y])[1])
    }
    
    #Define the vector of causal effects of intervention variables on var.y
    intervention.set.on.var.y <- sapply(intervention.set,adjusted.regression,y=var.y)
    
    #Define the required matrix of single intervention effects
    if (length(intervention.set)>1){
      intervention.set.on.intervention.set <- matrix(apply(expand.grid(intervention.set,intervention.set),
                                                           1,function(x) adjusted.regression(x[1],x[2])),nrow=length(intervention.set))
    }else{
      return(intervention.set.on.var.y)
    }

    joint.effect.fun  <- function(x){
      if(is.element(var.y,unlist(pasets[match(x,intervention.set)]))) return(0)
      x.temp  <- match(x,intervention.set)
      intervention.set.temp  <- intervention.set[-x.temp]
      #Intitialize the RR estimate as the single intervention effect of intervention.set on var.y
      RR.estimate <-  intervention.set.on.var.y[x.temp]
      #Define the vector of causal effects of intervention.set on other intervention variables
      x.temp.on.intervention.set.temp <- intervention.set.on.intervention.set[x.temp,-x.temp]
      #Define the vector of causal effects of "other" intervention variables on var.y
      intervention.set.on.var.y.temp <- intervention.set.on.var.y[-x.temp]
      #Define the required matrix of single intervention effects
      intervention.set.on.intervention.set.temp <- intervention.set.on.intervention.set[-x.temp,-x.temp]
      
      while (length(intervention.set.temp)>1){
        #update RR.estimate and the other things accounting for the elimination of the first entry of the current intervention set
        RR.estimate <- RR.estimate - x.temp.on.intervention.set.temp[1]*intervention.set.on.var.y.temp[1]
        intervention.set.on.var.y.temp <- intervention.set.on.var.y.temp[-1] - intervention.set.on.intervention.set.temp[-1,1] * intervention.set.on.var.y.temp[1]
        x.temp.on.intervention.set.temp <- x.temp.on.intervention.set.temp[-1] - x.temp.on.intervention.set.temp[1]*intervention.set.on.intervention.set.temp[1,-1]
        if (length(intervention.set.temp)>2)
          intervention.set.on.intervention.set.temp <- intervention.set.on.intervention.set.temp[-1,-1] - matrix(intervention.set.on.intervention.set.temp[-1,1],ncol=1)%*%matrix(intervention.set.on.intervention.set.temp[1,-1],nrow=1)
        intervention.set.temp <- intervention.set.temp[-1]
      }
      return(RR.estimate-x.temp.on.intervention.set.temp*intervention.set.on.var.y.temp)
    }
    
    return(sapply(intervention.set,joint.effect.fun))
    
}
