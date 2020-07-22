myrandomDAG2 <- function (n, fprob1, fprob2, ens, pwr=0, lB = 0.1, uB = 1, V = as.character(1:n)) 
{
  temp <- sum(((1:(n-1))^pwr)*(((n-1):1)^(1-pwr)))
  prob.global1 <- n*ens/(2*temp*fprob1)
  prob.global2 <- n*ens/(2*temp*fprob2)
  edL <- vector("list", n)
  nmbEdges <- 0L
  #panodes <- sample(1:(n-1),size=rbinom(1,n-1,fprob))
  for (i in seq_len(n - 2)) {
    if (i==1){
      temp <- rbinom(1,1,fprob1)
      prob <- prob.global1 * (i/(n-i))^pwr
    }
    else if(!is.element(i,unlist(edL[1:(i-1)]))){
      temp <- rbinom(1,1,fprob1)
      prob <- prob.global1 * (i/(n-i))^pwr
    }else{
      temp <- rbinom(1,1,fprob2)
      prob <- prob.global2 * (i/(n-i))^pwr
    }
    if (temp){
      listSize <- rbinom(1, n - i, prob)
      nmbEdges <- nmbEdges + listSize
      edgeList <- sample(seq(i + 1, n), size = listSize)
      weightList <- runif(length(edgeList), min = lB, max = uB)
      edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    else edL[[i]] <- list(edges = integer(0), weights = numeric(0))
  }
  if (nmbEdges > 0) {
    prob <- fprob1 * prob.global1 * (n-1)^pwr 
    edL[[n - 1]] <- if (rbinom(1, 1, prob) == 1) 
      list(edges = n, weights = runif(1, min = lB, max = uB))
    else list(edges = integer(0), weights = numeric(0))
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    names(edL) <- V
    new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  }
  else new("graphNEL", nodes = V, edgemode = "directed")
}