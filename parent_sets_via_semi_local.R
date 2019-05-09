parent.sets.via.semi.local <- function(x.pos,amat.cpdag){
  #amat.cpdag[which(amat.cpdag != 0)] <- 1
  amat.undir <- amat.cpdag*t(amat.cpdag) 
  amat.dir <- amat.cpdag - amat.undir 
  pasets.dir <- which(amat.dir[,x.pos]!=0)
  #g.undir <- as(amat.undir,"graphNEL")
  #conn.comp <- lapply(connectedComp(g.undir),as.integer)
  #conn.comp.imp <- conn.comp[sapply(conn.comp,function(comp) length(intersect(x.pos,comp))>0)]
  conn.comp.imp <- graph.dfs(graph=graph.adjacency(amat.undir,mode='undirected'),x.pos,unreachable=F)$order
  conn.comp.imp <- conn.comp.imp[!is.nan(conn.comp.imp)]
    nvar <- length(conn.comp.imp)
    if (nvar==1){
      all.pasets <- list(pasets.dir)
    }else{
      conn.comp.mat <- amat.undir[conn.comp.imp,conn.comp.imp]
      cat("conn.comp.nodes=",nvar,"conn.comp.edges=",sum(conn.comp.mat)/2,"\n")
      if (nvar>12) return(NULL)
      rownames(conn.comp.mat) <- conn.comp.imp
      x.temp <- match(x.pos,conn.comp.imp)
        rownames(conn.comp.mat) <- colnames(conn.comp.mat) <- 1:nvar
        all.extensions <- allDags(conn.comp.mat,conn.comp.mat,NULL)
        pa.fun <- function(r){
          amat.temp <- matrix(all.extensions[r,],nrow=nvar)
          return(c(conn.comp.imp[which(amat.temp[,x.temp]!=0)],pasets.dir))
        } 
        all.pasets <- lapply(1:nrow(all.extensions),pa.fun)
      
    }
  
  #all.pasets <- lapply(1:length(conn.comp.imp),extract.parent.sets.from.conn.comp)
  #idx <- expand.grid(lapply(1:length(all.pasets),function(i) 1:length(all.pasets[[i]])))
 # x.conn.comp <- unlist(lapply(1:length(conn.comp.imp),function(i) intersect(conn.comp.imp[[i]],x.pos)))
  #all.pasets <- lapply(1:nrow(idx),function(i) unlist(lapply(1:length(conn.comp.imp),
       #                                                      function(j) all.pasets[[j]][[idx[i,j]]]),recursive=F)[match(x.pos,x.conn.comp)])
  
  return(all.pasets)
}