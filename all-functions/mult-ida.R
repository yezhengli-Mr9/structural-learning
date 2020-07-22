multida <- function(x.pos,y.pos,mcov,graphEst=NULL,all.pasets=NULL,technique=c("RR","MCD")){
                                      
  if (is.null(all.pasets)){
    amat <- as(graphEst,"matrix")
    amat[which(amat != 0)] <- 1
    all.pasets <- extract.parent.sets(x.pos,amat)
    all.pasets <- all.pasets$pasets
  }else{
    if (length(x.pos)==1 & !is.list(all.pasets)) all.pasets <- list(all.pasets)
    if (length(x.pos)>1 & !is.list(all.pasets[[1]]))  all.pasets <- list(all.pasets)
  }
  if (length(y.pos)>1){
    return(lapply(y.pos,function(y) multida(x.pos,y,mcov,all.pasets=all.pasets,technique=technique)))
  }else{
    technique <- match.arg(technique)
    if (is.element(y.pos,x.pos)) return(matrix(0,nrow=length(x.pos),ncol=length(all.pasets)))
    joint.effects <- switch(technique,
                            RR = matrix(unlist(lapply(all.pasets,function(pasets) RR(mcov,x.pos,y.pos,pasets))),
                                        nrow=length(x.pos)),
                            MCD = matrix(unlist(lapply(all.pasets,function(pasets) MCD(mcov,x.pos,y.pos,pasets))),
                                         nrow=length(x.pos)))
    return(joint.effects)
  }
}