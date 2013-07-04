xform <- function(mat, scale=FALSE){
  mat0<-mat
  mat0[mat==0] <- NA
  smat <- scale(mat0, scale=scale)
  smat[is.na(smat)] <- 0
  smat
  }

makeMM <- function(FD, type, rescale=TRUE){
  if(type == 'zeroinf'){
    return(scale(exprs(FD)))
  } else if(type == 'hurdle'){
    return(cbind(xform(exprs(FD), scale=rescale), scale((exprs(FD)>0)*1, scale=rescale)))
  } else if(type=='dichotomous'){
    return( scale((exprs(FD[, freq(FD)<1])>0)*1, scale=rescale))
  } else if(type=='raw'){
    olayer <- layer(FD)
    layer(FD) <- 'lCount'
    model.mat <- scale(exprs(FD), scale=rescale)
    layer(FD) <- olayer
    return(model.mat)
  } else{
   stop('bad type') 
  }
}

N_COMPONENTS <- 30

## doPCA <- function(FD, type){
##     mm <- makeMM(FD, type)
##     fet0 <- SingleCellAssay::filter(FD, apply_filter=FALSE, filt_control=list(filter=FALSE))$fet0
##     prcomp(mm
    
    
## }
