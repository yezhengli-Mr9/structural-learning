myrandomDAG <- function (n, fprob, ens, lB = 0.1, uB = 1, V = as.character(1:n)) 
{
prob <- ens/((n-1)*fprob)
  edL <- vector("list", n)
  nmbEdges <- 0L
  panodes <- sample(1:(n-1),size=rbinom(1,n-1,fprob))
  for (i in seq_len(n - 2)) {
    if (is.element(i,panodes)){
      listSize <- rbinom(1, n - i, prob)
      nmbEdges <- nmbEdges + listSize
      edgeList <- sample(seq(i + 1, n), size = listSize)
      weightList <- runif(length(edgeList), min = lB, max = uB)
      edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    else edL[[i]] <- list(edges = integer(0), weights = numeric(0))
  }
  if (nmbEdges > 0) {
    edL[[n - 1]] <- if (is.element(n-1,panodes) & rbinom(1, 1, prob) == 1) 
      list(edges = n, weights = runif(1, min = lB, max = uB))
    else list(edges = integer(0), weights = numeric(0))
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    names(edL) <- V
    new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  }
  else new("graphNEL", nodes = V, edgemode = "directed")
}