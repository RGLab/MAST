xform <- function(mat, scale=FALSE){
  mat0<-mat
  mat0[mat==0] <- NA
  smat <- scale(mat0, scale=scale)
  smat[is.na(smat)] <- 0
  smat
  }
