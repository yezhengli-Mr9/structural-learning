parent.sets.via.local.IDA <- function(x,amat){
 # amat[which(amat != 0)] <- 1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  wgt.est <- (amat != 0)
  tmp <- wgt.est - t(wgt.est)
  tmp[which(tmp < 0)] <- 0
  wgt.unique <- tmp
  paset <- NULL
  pa1 <- which(wgt.unique[x, ] != 0)
  wgt.ambig <- wgt.est - wgt.unique
  pa2 <- which(wgt.ambig[x, ] != 0)
  if (length(pa2) == 0) {
    paset <- list(pa1)
  }
  else {
    pa2.f <- pa2
    pa2.t <- NA
    if (!has.new.coll(amat, amatSkel, x, pa1, pa2.t, 
                      pa2.f)) {
      paset <- list(pa1)
    }
    for (i2 in seq_along(pa2)) {
      pa2.f <- pa2[-i2]
      pa2.t <- pa2[i2]
      if (!has.new.coll(amat, amatSkel, x, pa1, pa2.t, pa2.f)) 
        paset <- c(paset,list(c(pa1,pa2.t)))
    }
    if (length(pa2) > 1) 
      for (i in 2:length(pa2)) {
        pa.tmp <- combn(pa2, i, simplify = TRUE)
        for (j in seq_len(ncol(pa.tmp))) {
          pa2.t <- pa.tmp[, j]
          pa2.f <- setdiff(pa2, pa2.t)
          if (!has.new.coll(amat, amatSkel, x, pa1, pa2.t, pa2.f))
            paset <- c(paset,list(c(pa1,pa2.t)))
        }
      }
  }
  
  return(paset)
}