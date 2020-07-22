extract.parent.sets.old <- function(x.pos,amat.cpdag){
  amat.cpdag[which(amat.cpdag != 0)] <- 1
  amat.undir <- amat.cpdag*t(amat.cpdag) 
  amat.dir <- amat.cpdag - amat.undir 
  pasets.dir <- lapply(x.pos,function(x) which(amat.dir[,x]!=0))
  g.undir <- as(amat.undir,"graphNEL")
  conn.comp <- lapply(connectedComp(g.undir),as.integer)
  conn.comp.imp <- conn.comp[sapply(conn.comp,function(comp) length(intersect(x.pos,comp))>0)]
  chrordality.test.fun <- function(comp) is.chordal(graph.adjacency(amat.undir[comp,comp],mode="undirected"),fillin=F)$chordal
  chordal <- sapply(conn.comp.imp, chrordality.test.fun)

all.locally.valid.parents.undir <- function(amat,x){ #x must be a scaler
  pa.dir <- pasets.dir[[which(x.pos==as.integer(rownames(amat)[x]))]]
  paset <- list(pa.dir)
  pa <- which(amat[,x] != 0) #cannot be a null set
  if (length(pa)==1){
    paset <- c(paset,list(c(as.integer(rownames(amat)[pa]),pa.dir))) 
  }else{
    for (i in 1:length(pa)){
      pa.tmp <- combn(pa, i, simplify = TRUE)
      n.comb <- ncol(pa.tmp)
      for (j in 1:n.comb) {
        pa.t <- pa.tmp[, j]
        new.coll <- FALSE
        if (length(pa.t)>1){
          tmp <- amat[pa.t,pa.t]
          diag(tmp) <- 1
          new.coll <- (min(tmp)==0)
        }
        if (!new.coll) paset <- c(paset,list(c(as.integer(rownames(amat)[pa.t]),pa.dir)))
      }
    }
  }
  
  return(paset)
  }

  extract.parent.sets.from.conn.comp <- function(i){
    all.nodes <- conn.comp.imp[[i]]
    nvar <- length(all.nodes)
    if (nvar==1){
      pasets.comp <- list(pasets.dir[match(all.nodes,x.pos)])
    }else{
      conn.comp.mat <- amat.undir[all.nodes,all.nodes]
      rownames(conn.comp.mat) <- all.nodes
      x.temp <- intersect(all.nodes,x.pos)
      x.temp2 <- match(x.temp,all.nodes)
      if(chordal[i]){
        rownames(conn.comp.mat) <- colnames(conn.comp.mat) <- 1:nvar
        all.extensions <- allDags(conn.comp.mat,conn.comp.mat,NULL)
        pa.fun <- function(amat,j) c(all.nodes[which(amat[,x.temp2[j]]!=0)],pasets.dir[[match(x.temp[j],x.pos)]])
        parent.sets.fun <- function(r) lapply(1:length(x.temp),pa.fun,amat=matrix(all.extensions[r,],nrow=nvar))
        pasets.comp <- lapply(1:nrow(all.extensions),parent.sets.fun)
      }else{
        pasets.comp <- lapply(x.temp2,all.locally.valid.parents.undir,amat=conn.comp.mat)
        idx <- expand.grid(lapply(1:length(x.temp),function(j) 1:length(pasets.comp[[j]])))
        pasets.comp <- lapply(1:nrow(idx),function(r) lapply(1:length(x.temp), function(j) pasets.comp[[j]][[idx[r,j]]]))
      }
      
    }
    return(pasets.comp)
  }

all.pasets <- lapply(1:length(conn.comp.imp),extract.parent.sets.from.conn.comp)
idx <- expand.grid(lapply(1:length(all.pasets),function(i) 1:length(all.pasets[[i]])))
x.conn.comp <- unlist(lapply(1:length(conn.comp.imp),function(i) intersect(conn.comp.imp[[i]],x.pos)))
all.pasets <- lapply(1:nrow(idx),function(i) unlist(lapply(1:length(conn.comp.imp),
                                                function(j) all.pasets[[j]][[idx[i,j]]]),recursive=F)[match(x.pos,x.conn.comp)])

return(list(pasets=all.pasets,n.nonchordal=sum(!chordal)))
}