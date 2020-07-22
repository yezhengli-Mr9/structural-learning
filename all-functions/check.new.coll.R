check.new.coll  <- function (amat, amatSkel, x, pa1, pa2.t, pa2.f){
  res <- FALSE
  if ((length(pa2.t) > 0) && any(!is.na(pa2.t))) {
    if ((length(pa1) > 0) && any(!is.na(pa1))) {
      res <- (min(amatSkel[pa1, pa2.t]) == 0)
    }
    if ((length(pa2.t) > 1) && (!res)) {
      tmp <- amatSkel[pa2.t, pa2.t]
      diag(tmp) <- 1
      res2 <- (min(tmp) == 0)
      res <- (res | res2)
    }
  }
  if (!res && ((length(pa2.f) > 0) && any(!is.na(pa2.f)))) {
    amatTmp <- amat
    amatTmp <- amatTmp - t(amatTmp)
    amatTmp[amatTmp < 0] <- 0
    tmp <- amatTmp[pa2.f, , drop = FALSE]
    papa <- setdiff(which(apply(tmp, 2, sum) != 0), x)
    if (length(papa) > 0) 
      res <- res | (min(amatSkel[x, papa]) == 0)
  }
  return(res)
}