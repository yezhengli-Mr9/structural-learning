#source('http://bioconductor.org/biocLite.R')
#biocLite('Rgraphviz')
library(Rgraphviz)
library(pcalg)
#source('http://bioconductor.org/biocLite.R')
#biocLite('RBGL')
library(RBGL)
library(mvtnorm)
library(igraph)
#biocLite('multicore')
library(multicore)
#library(symmoments)
library(corpcor)
library(bdsmatrix)
library(Matrix)
library(pracma)
library(huge)
source("myrandomDAG.R")
source("nonGaussian_simulation.R")
source("run-simulation.R")
#Set working directory
setwd("/Users/yezheng/github/structural-learning")
#Load all required functions
all.functions <- sapply(list.files(pattern="[.]R$", path="all-functions", full.names=TRUE), source)
seed <- 123
p <- 1000
ens <- 4
fprob <- 1
nsample <- 1000
alpha <- 0.01
fprob <- 1/5
nullparents <- F
for (fprob in c(1,1/5)){
  for (nullparents in c(T,F)){
    predictions <- function(j){
      cat("\n",j)
      set.seed(seed*j)
      #g <- randomDAG(p,ens/(p-1),lB=0.2,uB=1.2)
      g <- myrandomDAG(p,fprob=fprob,ens=ens,lB=0.2,uB=1.2)
      print(g)
      amat <- as(g,"matrix")
      #   amat <- as(huge.generator(d=p,graph="band",g=2)$theta,'matrix')*matrix(runif(p*p,min=0.2,max=1.2),nrow=p)
      #   g <- as(amat,'graphNEL')
      #   print(g)
      #   idx.perm <- randperm(1:p)
      #   #amat <- amat + t(amat)
      #   amat <- amat[idx.perm,idx.perm]
      #   amat[lower.tri(amat)] <- 0
      #   rownames(amat) <- colnames(amat) <- as.character(1:p)
      #   g <- as(amat,"graphNEL")
      #g.cpdag <- dag2cpdag(g)
      #amat.cpdag <- as(g.cpdag,"matrix")
      error.cov.mat <- diag(p)
      true.cov.mat<- solve(diag(p)-t(amat))%*%error.cov.mat%*%solve(diag(p)-amat)
      #   nparents <- unlist(lapply(1:p,function(x) length(which(amat[,x]!=0))))
      #   imp.nodes <- which(nparents>0)
      #   all.intervention.set <- combn(imp.nodes,2)
      #   all.intervention.set <- all.intervention.set[,sample(1:ncol(all.intervention.set),10)]
      cor.pairs <- which(abs(true.cov.mat)>0.001,arr.ind=T)
      nparents <- unlist(lapply(1:p,function(x) length(which(amat[,x]!=0))))
      nchildren <- unlist(lapply(1:p,function(x) length(which(amat[x,]!=0))))
      idx <- which(nparents>0)
      idx <- unique(c(idx,which(nchildren>0)))
      all.intervention.set <- combn(idx,2)
      all.intervention.set <- all.intervention.set[,sample(1:ncol(all.intervention.set),10)]
      #all.intervention.set <- unique(cor.pairs[,1])
      #generate data
      idx <- sample(c(T,F),size=p,replace=T)
      mix.error.fun <- function(i) if(idx[i]) rt(nsample,df=3)/sqrt(3) else rnorm(nsample)
      data <- sapply(1:p,mix.error.fun)
      data <- data%*%solve(diag(p)-amat)
      data <- apply(data,2,function(x) x -mean(x))
      samp.cov.mat <- cov(data)
      #obtain graph estimate
      suffStat.data <- list(C=cor(data), n=nsample)
      #suffStat.data <- list(C=cor(data), n=nsample)
      indepTest.data <- gaussCItest
      cat("Running PC...")
      if (!nullparents){
        graphEst <- pc(suffStat.data, indepTest.data, p=p, alpha=alpha,u2pd="relaxed")@graph
      }
      
      cat("  Done.\n")
      #all.intervention.set <- combn(p,2)[,sample(1:(p*(p-1)/2),20)]
      #all.intervention.set <- combn(p,2)
      
      joint.effects.given.intervention.set <- function(i){
        #cat("\n",i)
        intervention.set <- all.intervention.set[,i]
        all.y <- union(cor.pairs[which(cor.pairs[,1]==intervention.set[1]),2],
                       cor.pairs[which(cor.pairs[,1]==intervention.set[2]),2])
        #all.y <- cor.pairs[which(cor.pairs[,1]==intervention.set),2]
        pasets <- lapply(intervention.set,function(x) which(amat[,x]!=0))
        #all.y <- (1:p)[-intervention.set]
        all.y <- setdiff(all.y,intervention.set)
        true.joint.effects.DAG <- as.vector(do.call("cbind",lapply(all.y,
                                                                   function(var.y) RRC(true.cov.mat,intervention.set,var.y,pasets))))
        
        #all.pasets.true <- extract.parent.sets(intervention.set,amat.cpdag,isCPDAG=TRUE)
        all.pasets.true <- list(list(NULL,NULL))
        true.joint.effects.CPDAG <- lapply(all.y,function(var.y) 
          jointIda(intervention.set,var.y,true.cov.mat,all.pasets=all.pasets.true,technique="RRC"))
        true.joint.effects.MinAbs <- as.vector(do.call("cbind",lapply(true.joint.effects.CPDAG,function(x) apply(abs(x),1,min))))
        true.joint.effects.Aver <- as.vector(do.call("cbind",lapply(true.joint.effects.CPDAG,function(x) apply(x,1,mean))))
        true.joint.effects <- cbind(true.joint.effects.DAG,true.joint.effects.MinAbs,true.joint.effects.Aver)
        colnames(true.joint.effects) <- c("DAG","MinAbs","Aver")
        
        
        all.pasets.estimated <- if(nullparents) list(list(NULL,NULL)) else extract.parent.sets(intervention.set,as(graphEst,"matrix"))
        #all.pasets.estimated <- list(list(NULL,NULL))
        estimated.joint.effects.IPW <- lapply(all.y,function(var.y) 
          matrix(unlist(lapply(all.pasets.estimated,function(pasets) IPW(data,intervention.set,var.y,pasets))),nrow=length(intervention.set)))
        estimated.joint.effects.MinAbs.IPW <- as.vector(do.call("cbind",lapply(estimated.joint.effects.IPW,function(x) apply(abs(x),1,min))))
        estimated.joint.effects.Aver.IPW <- as.vector(do.call("cbind",lapply(estimated.joint.effects.IPW,function(x) apply(x,1,mean))))
        
        estimated.joint.effects.MCD <- lapply(all.y,function(var.y) 
          jointIda(intervention.set,var.y,samp.cov.mat,all.pasets=all.pasets.estimated,technique="MCD"))
        estimated.joint.effects.MinAbs.MCD <- as.vector(do.call("cbind",lapply(estimated.joint.effects.MCD,function(x) apply(abs(x),1,min))))
        estimated.joint.effects.Aver.MCD <- as.vector(do.call("cbind",lapply(estimated.joint.effects.MCD,function(x) apply(x,1,mean))))
        
        estimated.joint.effects.RRC <- lapply(all.y,function(var.y) 
          jointIda(intervention.set,var.y,samp.cov.mat,all.pasets=all.pasets.estimated,technique="RRC"))
        estimated.joint.effects.MinAbs.RRC <- as.vector(do.call("cbind",lapply(estimated.joint.effects.RRC,function(x) apply(abs(x),1,min))))
        estimated.joint.effects.Aver.RRC <- as.vector(do.call("cbind",lapply(estimated.joint.effects.RRC,function(x) apply(x,1,mean))))
        
        estimated.joint.effects <- list(MinAbs=NULL,Aver=NULL)
        estimated.joint.effects$MinAbs <- cbind(estimated.joint.effects.MinAbs.IPW,estimated.joint.effects.MinAbs.MCD,estimated.joint.effects.MinAbs.RRC)
        colnames(estimated.joint.effects$MinAbs) <- c("IPW","MCD","RRC")
        estimated.joint.effects$Aver <- cbind(estimated.joint.effects.Aver.IPW,estimated.joint.effects.Aver.MCD,estimated.joint.effects.Aver.RRC)
        colnames(estimated.joint.effects$Aver) <- c("IPW","MCD","RRC")
        
        return(list(true.joint.effects=true.joint.effects,estimated.joint.effects=estimated.joint.effects))
      }
      
      
      result <- mclapply(1:ncol(all.intervention.set),joint.effects.given.intervention.set,mc.cores=4)
      result <- do.call("rbind",result)
      true.joint.effects <- do.call("rbind",result[,1])
      #true.joint.effects <- unlist(result[,1])
      estimated.joint.effects <- do.call("rbind",result[,2])
      estimated.joint.effects.MinAbs <- do.call("rbind",estimated.joint.effects[,1])
      estimated.joint.effects.Aver <- do.call("rbind",estimated.joint.effects[,2])
      
      return(list(goldstd.effects=abs(true.joint.effects),
                  estimates.Aver=list(IPW=estimated.joint.effects.Aver[,1],
                                      MCD=estimated.joint.effects.Aver[,2],
                                      RRC=estimated.joint.effects.Aver[,3]),
                  estimates.MinAbs=list(IPW=estimated.joint.effects.MinAbs[,1],
                                        MCD=estimated.joint.effects.MinAbs[,2],
                                        RRC=estimated.joint.effects.MinAbs[,3])))
    }
    
    
    
    all.predictions <- lapply(1:5,predictions)
    
    #  filename <- paste("Results/estimated-CPDAG-p=",p,"-ens=",ens,"-nsample=",nsample,"-alpha=",alpha,"-nonGaussian.Rdata",sep="")
    #  save(file=filename,all.predictions)
    
    #load("Results/estimated-CPDAG-p=1000-ens=4-nsample=1000-alpha=0.05-nonGaussian.Rdata")
    
    set.seed(10)
    # for (i in 1:5){
    #   #idx <- 1:20
    #   #idx <- unlist(lapply(1:10, function(x) (x-1)*396 + 1:396))
    #   idx <- unlist(lapply(1:10, function(x) (x-1)*996 + 2*sample(1:498,100) - 1))
    #   idx <- as.vector(rbind(idx,idx+1))
    #   all.predictions[[i]]$goldstd.effects <- all.predictions[[i]]$goldstd.effects[idx]
    #   all.predictions[[i]]$estimates.Aver$IPW <- all.predictions[[i]]$estimates.Aver$IPW[idx]
    #   all.predictions[[i]]$estimates.Aver$MCD <- all.predictions[[i]]$estimates.Aver$MCD[idx]
    #   all.predictions[[i]]$estimates.Aver$RRC <- all.predictions[[i]]$estimates.Aver$RRC[idx]
    # }
    
    DK.effects.goldstd <- lapply(all.predictions,function(x) as.vector(x$goldstd.effects[,1]))
    
    #predictions.goldstd <- lapply(all.predictions,function(x) as.vector(x$goldstd.effects[,3]))
    
    predictions.IPW <- lapply(all.predictions,function(x) order(abs(as.vector(x$estimates.Aver$IPW)),decreasing=T))
    predictions.MCD <- lapply(all.predictions,function(x) order(abs(as.vector(x$estimates.Aver$MCD)),decreasing=T))
    predictions.RR <- lapply(all.predictions,function(x) order(abs(as.vector(x$estimates.Aver$RRC)),decreasing=T))
    
    # predictions.IPW <- lapply(all.predictions,function(x) as.vector(x$estimates.MinAbs$IPW))
    # predictions.MCD <- lapply(all.predictions,function(x) as.vector(x$estimates.MinAbs$MCD))
    # predictions.RR <- lapply(all.predictions,function(x) as.vector(x$estimates.MinAbs$RRC))
    
    neffects <- unlist(lapply(DK.effects.goldstd,function(x) length(x)))
    #target <- ceiling(0.01*neffects)
    #target <- 100
    target <- lapply(DK.effects.goldstd,function(x) ceiling(length(which(x>0.0001))))
    #target <- lapply(1:5,function(i) min(ceiling(length(which(DK.effects.goldstd[[i]]>0.0001))),
    #                                  ceiling(0.1*neffects[[i]])))
    target.set <- lapply(1:5,function(i) order(abs(DK.effects.goldstd[[i]]),decreasing=T)[1:target[[i]]])
    
    set.seed(123)
    random.orderings <- lapply(1:5,function(i) sapply(1:5,function(x) randperm(neffects[i],neffects[i])))
    
    get.roc <- function(DK.effects.est,target.set){
      # order.est <- order(abs(DK.effects.est),decreasing=T)
      order.est <- DK.effects.est
      true.positive <- sapply(1:length(DK.effects.est),function(x) length(intersect(target.set,order.est[1:x])))
      false.positive <- 1:length(DK.effects.est) - true.positive
      true.positive <- true.positive/length(target.set)
      false.positive <- false.positive/(length(DK.effects.est)-length(target.set))
      return(cbind(false.positive,true.positive))
    }
    
    get.pr <- function(DK.effects.est,target.set){
      # order.est <- order(abs(DK.effects.est),decreasing=T)
      order.est <- DK.effects.est
      true.positive <- sapply(1:length(DK.effects.est),function(x) length(intersect(target.set,order.est[1:x])))
      precision <- true.positive/(1:length(DK.effects.est))
      recall <- true.positive/length(target.set)
      return(cbind(recall,precision))
    }
    
    all.ROC.RR <- do.call('rbind',mclapply(1:5,function(i) get.roc(predictions.RR[[i]],target.set[[i]]),mc.cores=3))
    all.ROC.MCD <- do.call('rbind',mclapply(1:5,function(i) get.roc(predictions.MCD[[i]],target.set[[i]]),mc.cores=3))
    all.ROC.IPW <- do.call('rbind',mclapply(1:5,function(i) get.roc(predictions.IPW[[i]],target.set[[i]]),mc.cores=3))
    #all.ROC.goldstd <- do.call('rbind',mclapply(1:5,function(i) get.roc(predictions.goldstd[[i]],target.set[[i]]),mc.cores=3))
    all.ROC.random <- do.call('rbind',mclapply(1:5, function(j) do.call('rbind',
                          lapply(1:5,function(i) get.roc(random.orderings[[i]][,j],target.set[[i]]))),mc.cores=3))
    
    #all.ROC <- rbind(all.ROC.IPW,all.ROC.MCD,all.ROC.RR,all.ROC.goldstd,all.ROC.random)
    all.ROC <- rbind(all.ROC.IPW,all.ROC.MCD,all.ROC.RR,all.ROC.random)
    
    all.PR.RR <- do.call('rbind',mclapply(1:5,function(i) get.pr(predictions.RR[[i]],target.set[[i]]),mc.cores=3))
    all.PR.MCD <- do.call('rbind',mclapply(1:5,function(i) get.pr(predictions.MCD[[i]],target.set[[i]]),mc.cores=3))
    all.PR.IPW <- do.call('rbind',mclapply(1:5,function(i) get.pr(predictions.IPW[[i]],target.set[[i]]),mc.cores=3))
    #all.PR.goldstd <- do.call('rbind',mclapply(1:5,function(i) get.pr(predictions.goldstd[[i]],target.set[[i]]),mc.cores=3))
    all.PR.random <- do.call('rbind',mclapply(1:5, function(j) do.call('rbind',
           lapply(1:5,function(i) get.pr(random.orderings[[i]][,j],target.set[[i]]))),mc.cores=3))
    
    #all.PR <- rbind(all.PR.IPW,all.PR.MCD,all.PR.RR,all.PR.goldstd,all.PR.random)
    all.PR <- rbind(all.PR.IPW,all.PR.MCD,all.PR.RR,all.PR.random)
    all.values <- round(rbind(all.ROC,all.PR),5)
    
    Network.id <- unlist(lapply(1:5,function(i) paste("Network",i,sep=" ")))
    methods <- c("IPW","MCD","RRC",unlist(lapply(1:5,function(i) paste("Random.",i,sep=""))))
    
    all.values <- data.frame(values.x =all.values[,1], values.y= all.values[,2], Network.id= factor(rep(rep(Network.id,times=neffects),8*2)), 
                             curve.type = factor(rep(c("ROC","Precision-Recall"),each=sum(neffects)*8),
                                                 levels=c("Precision-Recall","ROC") ), Method = rep(rep(methods,each=sum(neffects)),2) )
    library(lattice)
    library(latticeExtra)
    
    temp1 <- if (fprob==0.2) "hub" else "ER"
    temp2 <- if (nullparents) "null_pasets" else "estimated_pasets"
    
    filename <- paste("ROC_PR_",temp1,"_",temp2,"_PC.pdf",sep="")
    #postscript("whatever.eps")
    pdf(filename,width=7.8,height=3.9,pointsize=6,compress=T)
    par(mfrow = c(2,5),
        oma = 0.1*c(1,1,1,1),
        mar = 0.1*c(1,1,1,1))
    useOuterStrips(xyplot(values.y~values.x|Network.id * curve.type, data= all.values, groups = Method,xlab='',ylab='',asp=1,
                          panel = panel.superpose,lty=c(4,5,1,1,1,1,0,0),col=c('black','red','darkblue',rep('gray',5)),type='s',
                          scales=list(tick.number=2,at=seq(0.1,1,0.2),lables=seq(0.1,1,0.2),cex=0.7),
                          key=list(x=.109,y=.475,corner=c(0,0),border=T,
                                   lines=list(lty=c(4,5,1,1),col=c('black','red','darkblue','gray')),
                                   text=list(c("IPW","MCD","RRC","Random")),cex=0.4,size=1.5),
                          par.settings = list(strip.background=list(col="lightgrey"))))
    dev.off()
     
    
  }
}
