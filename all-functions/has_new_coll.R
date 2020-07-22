has.new.coll <- function(amat,amatSkel, x, pa1, pa2.t, pa2.f) {
  ## Check if undirected edges that are pointed to x create a new v-structure
  ## Additionally check, if edges that are pointed away from x create
  ## new v-structure; i.e. x -> pa <- papa would be problematic
  ## pa1 are definit parents of x
  ## pa2 are undirected "parents" of x
  ## pa2.t are the nodes in pa2 that are directed towards pa2
  ## pa2.f are the nodes in pa2 that are directed away from pa2
  ## Value is TRUE, if new collider is introduced
  res <- FALSE
  if (length(pa2.t) > 0 && !all(is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if no, there is a new collider
    if (length(pa1) > 0 && !all(is.na(pa1))) {
      res <- min(amatSkel[pa1, pa2.t]) == 0 ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if (!res && length(pa2.t) > 1) {
      A2 <- amatSkel[pa2.t,pa2.t]
      diag(A2) <- 1
      res <- min(A2) == 0 ## TRUE if new collider
    }
  }
  if (!res && length(pa2.f) > 0 && !all(is.na(pa2.f))) {
    ## consider here only the DIRECTED Parents of pa2.f
    ## remove undirected edges
    A <- amat-t(amat)
    A[A<0] <- 0
    ## find parents of pa2.f
    cA <- colSums(A[pa2.f,,drop=FALSE])
    papa <- setdiff(which(cA != 0), x)
    ## if any node in papa is not directly connected to x, there is a new
    ## collider
    if (length(papa) > 0)
      res <- min(amatSkel[x,papa]) == 0 ## TRUE if new collider
  }
  res
}