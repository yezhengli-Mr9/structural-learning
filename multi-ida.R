multi.ida <- function(x.pos,y.pos,mcov,graphEst,method=c("RR","MCD","Both"),summary.measure=c("None","MinAbs","Aver","Both")){
  method <- match.arg(method)
  summary.measure <- match.arg(summary.measure)
  joint.effects <- list(RR=NULL,MCD=NULL)
  if (is.element(y.pos,x.pos)){
    joint.effects$RR <- if (method == "RR" | method == "Both") rep(0,length(x.pos))
    joint.effects$MCD <- if (method == "MCD" | method == "Both") rep(0,length(x.pos))
    return(joint.effects)
  }
  amat <- as(graphEst,"matrix")
  amat[which(amat != 0)] <- 1
  all.pasets <- extract.parent.sets(x.pos,amat)
  all.pasets <- all.pasets$pasets
#   imp.var <- unique(c(x.pos,unlist(all.pasets),y.pos))
#   if (length(imp.var) < nrow(mcov)){
#     mcov <- mcov[imp.var,imp.var]
#     x.pos <- match(x.pos,imp.var)
#     y.pos <- match(y.pos,imp.var)
#     pasets <- lapply(pasets,function(x) match(unlist(x),imp.var))
#   }

  if (method == "RR" | method == "Both"){
    joint.effects$RR <- matrix(unlist(lapply(all.pasets,function(pasets) RR(mcov,x.pos,y.pos,pasets))),
                                  nrow=length(x.pos))
    if (!(summary.measure=="None"))
      joint.effects$RR <- switch(summary.measure,
                                 MinAbs = apply(abs(joint.effects$RR),1,min),
                                 Aver = apply(abs(joint.effects$RR),1,mean),
                                 Both = list(MinAbs=apply(abs(joint.effects$RR),1,min), Aver=apply(joint.effects$RR,1,mean)))
   }
    
   if (method == "MCD" | method == "Both"){
     joint.effects$MCD <- matrix(unlist(lapply(all.pasets,function(pasets) MCD(mcov,x.pos,y.pos,pasets))),
                                nrow=length(x.pos))
     if (!(summary.measure=="None"))
       joint.effects$MCD <- switch(summary.measure,
                                  MinAbs = apply(abs(joint.effects$MCD),1,min),
                                  Aver = apply(abs(joint.effects$MCD),1,mean),
                                  Both = list(MinAbs=apply(abs(joint.effects$MCD),1,min), Aver=apply(joint.effects$MCD,1,mean)))
   }
return(joint.effects)
}