IPW <- function(data,intervention.set,var.y,pasets){
  
  if (is.element(var.y,intervention.set)) return(rep(0,length(intervention.set)))
  data <- apply(data,2,function(x) x -mean(x))
  weights.fun <- function(x){
    pa <- unlist(pasets[match(x,intervention.set)])
    if (length(pa)==0)
      return(matrix(1,nrow=nrow(data),ncol=1))
    fit <- lm(data[,x] ~ 0+data[,pa])
    return(exp(dnorm(data[,x],sd=sd(data[,x]),log=TRUE) - 
                 dnorm(data[,x],mean= data[,pa]%*%as.matrix(fit$coeff[1:length(pa)],nrow=length(pa)),
                       sd=sd(fit$resid),log=TRUE)))
  }
  
  weights <- apply(sapply(matrix(intervention.set,nrow=1),weights.fun),1,prod)

ipw.estimate <- lm(data[,var.y] ~ 0+data[,intervention.set], weights=weights)$coeff
names(ipw.estimate) <- NULL
for (i in 1:length(intervention.set))
  if(is.element(var.y,unlist(pasets[match(intervention.set[i],intervention.set)]))) ipw.estimate[i] <- 0
return(ipw.estimate)
}