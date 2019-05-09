jointIda <- function(x.pos,y.pos,mcov,graphEst=NULL,all.pasets=NULL,technique=c("RRC","MCD")){
  
  if (is.null(all.pasets)){
    amat <- as(graphEst,"matrix")
    amat[which(amat != 0)] <- 1
    all.pasets <- extract.parent.sets(x.pos,amat)
  }else{
    correct.format <- TRUE
    if (!is.list(all.pasets)) correct.format <- FALSE
    for (i in 1:length(all.pasets)){
      if (length(all.pasets[[i]]) != length(x.pos)) correct.format <- FALSE
    }
      if (!correct.format) stop("all.pasets is not given in an appropriate format.")
  }
  if (length(y.pos)>1){
    return(lapply(y.pos,function(y) jointIda(x.pos,y,mcov,all.pasets=all.pasets,technique=technique)))
  }else{
    technique <- match.arg(technique)
    if (is.element(y.pos,x.pos)) return(matrix(0,nrow=length(x.pos),ncol=length(all.pasets)))
    joint.effects <- switch(technique,
                            RRC = matrix(unlist(lapply(all.pasets,function(pasets) RRC(mcov,x.pos,y.pos,pasets))),
                                        nrow=length(x.pos)),
                            MCD = matrix(unlist(lapply(all.pasets,function(pasets) MCD(mcov,x.pos,y.pos,pasets))),
                                         nrow=length(x.pos)))
    return(joint.effects)
  }
}