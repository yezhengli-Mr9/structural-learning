run.simulation <- function(simulation.parameters,ncores=1,known.DAG=F,known.CPDAG=F){
  block.size <- simulation.parameters$block.size
  nblocks <- simulation.parameters$nblocks
  p <- nblocks*block.size   #Number of nodes
  ens <- simulation.parameters$ens #Expected neighborhood size
  nintervention <- simulation.parameters$nintervention #number of joint interventions
  nsample <- simulation.parameters$nsample #sample size for each iteration
  niter <- simulation.parameters$niter # number of iterations
  
  cat("\n Running simulation with\n Number of variables =",p,"\n Expected neighborhood size =",ens,"\n and sample size =",nsample,"...")
  
  run.iteration <- function(i){
    set.seed(simulation.parameters$seed*i)
    nrep <- 1
    repeat{
      amat <- as.matrix(bdiag(lapply(1:nblocks,function(i) as(randomDAG(block.size,ens/(block.size-1),lB=0.2,uB=1.2),"matrix"))))
      rownames(amat) <- as.character(1:p); colnames(amat) <- as.character(1:p)
      #error.cov.mat <- diag(1+runif(p))
      error.cov.mat <- diag(p)
      true.cov.mat<- solve(diag(p)-t(amat))%*%error.cov.mat%*%solve(diag(p)-amat)
      #true.cov.mat <- cov2cor(true.cov.mat)
      g<- as(amat!=0,"graphNEL")
      
      #Choose intervention nodes and var.y
      block <- (sample.int(n=nblocks,size=1)-1)*block.size + (1:block.size)
      intervention.set <-  sample(x=block,size=nintervention,replace=FALSE)
      pasets <- lapply(intervention.set,function(x) which(amat[,x]!=0))
      imp.var <- unique(c(intervention.set,unlist(pasets)))
      if (length(imp.var)==block.size-1){
        var.y <- setdiff(block,imp.var)
      }else{
        var.y <- if(length(imp.var)<block.size) sample(x=setdiff(block,imp.var),size=1) else sample(x=setdiff(block,intervention.set),size=1)
      }
      nrep <- nrep+1
      if (length(imp.var)<block.size | nrep > 100) break
    }
    true.joint.effects.DAG <- RRC(true.cov.mat,intervention.set,var.y,pasets)
    
    if (known.DAG){
      true.causal.effect <- RRC(true.cov.mat,intervention.set,var.y,pasets)
    }else{
      g.cpdag <- dag2cpdag(g)
      #true.joint.effects <- multi.ida(intervention.set,var.y,true.cov.mat,g.cpdag,"RRC","Both")$RRC
      true.joint.effects <- jointIda(intervention.set,var.y,true.cov.mat,graphEst=g.cpdag,technique="RRC")
      true.joint.effects<- list(MinAbs=apply(abs(true.joint.effects),1,min),
                                          Aver=apply(true.joint.effects,1,mean))
    }
    #Generate data
    #data <- rmvnorm(nsample,mean=rep(0,p),sigma=true.cov.mat)
    idx <- sample(c(T,F),size=p,replace=T)
    mix.error.fun <- function(i) if(idx[i]) rt(nsample,df=3)/sqrt(3) else rnorm(nsample)
    data <- sapply(1:p,mix.error.fun)
    data <- data%*%solve(diag(p)-amat)
    data <- apply(data,2,function(x) x -mean(x))
    samp.cov.mat <- cov(data)
    estimated.joint.effects <- list(IPW=NULL,RRC=NULL,MCD=NULL)
    if (known.DAG){
      estimated.joint.effects$IPW  <- IPW(data,intervention.set,var.y,pasets)
      estimated.joint.effects$RRC  <- RRC(samp.cov.mat,intervention.set,var.y,pasets)
      estimated.joint.effects$MCD  <- MCD(samp.cov.mat,intervention.set,var.y,pasets)
    }else{
      if (known.CPDAG){
        graphEst <- g.cpdag
      }else{
        #suffStat.data <- list(C=cor(data), n=nsample)
        suffStat.data <- list(C=2*sin(cor(data,method="spearman")*pi/6), n=nsample)
        indepTest.data <- gaussCItest
        graphEst <- pc(suffStat.data, indepTest.data, p=p, alpha=simulation.parameters$alpha,u2pd="relaxed")@graph
      }
      all.pasets <- extract.parent.sets(intervention.set,as(graphEst,"matrix"))
      #all.chordal <- if (!known.CPDAG) (sum(all.pasets$n.nonchordal) == 0)
     # all.pasets <- all.pasets$pasets
      estimated.joint.effects$IPW <- matrix(unlist(lapply(all.pasets,function(pasets) IPW(data,intervention.set,var.y,pasets))),
                                                                     nrow=length(intervention.set))
      #temp <- multi.ida(intervention.set,var.y,samp.cov.mat,graphEst,"Both","Both")
      estimated.joint.effects$RRC  <- jointIda(intervention.set,var.y,samp.cov.mat,all.pasets=all.pasets,technique="RRC")
      estimated.joint.effects$MCD  <- jointIda(intervention.set,var.y,samp.cov.mat,all.pasets=all.pasets,technique="MCD")
      estimated.joint.effects$IPW <- list(MinAbs=apply(abs(estimated.joint.effects$IPW),1,min),
                                          Aver=apply(estimated.joint.effects$IPW,1,mean))
      estimated.joint.effects$RRC <- list(MinAbs=apply(abs(estimated.joint.effects$RRC),1,min),
                                          Aver=apply(estimated.joint.effects$RRC,1,mean))
      estimated.joint.effects$MCD <- list(MinAbs=apply(abs(estimated.joint.effects$MCD),1,min),
                                          Aver=apply(estimated.joint.effects$MCD,1,mean))
    }
    output <- list(true.joint.effects.DAG=true.joint.effects.DAG,true.joint.effects=true.joint.effects,estimated.joint.effects=estimated.joint.effects)
    #if (!known.CPDAG) output <- append(output,list(all.chordal=all.chordal))
    return(output)
  }

# for (i in 1:50){
#   print(i)
#   a <- run.iteration(i)
# }
  
  result.all <- mclapply(1:niter,run.iteration,mc.cores=ncores)
  result.all  <- do.call(rbind,result.all)
  true.joint.effects.DAG <- do.call(rbind,result.all[,1])
  true.joint.effects <- do.call(rbind,result.all[,2])
  estimated.joint.effects <- do.call(rbind,result.all[,3])
  estimated.joint.effects <- list(IPW=do.call(rbind,estimated.joint.effects[,1]),
                                  RRC=do.call(rbind,estimated.joint.effects[,2]),MCD=do.call(rbind,estimated.joint.effects[,3]))
  output <- list(true.joint.effects.DAG=true.joint.effects.DAG,true.joint.effects=true.joint.effects,estimated.joint.effects=estimated.joint.effects)
  #if (!known.CPDAG) output <- append(output,list(all.chordal=unlist(result.all[,3])))
  return(output)
  
}
