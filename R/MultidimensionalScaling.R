##' Make matrix of continuous expression values, orthogonal to discrete
##'
##' This centers each column of \code{mat} around the mean of its non-zero values.
##' @param mat matrix (such as produced by exprs)
##' @param scale should the columns also be scaled to have unit variance
##' @export
xform <- function(mat, scale=FALSE){
    mat0<-mat
    mat0[mat==0] <- NA
    smat <- scale(mat0, scale=scale)
    smat[is.na(smat)] <- 0
    smat
}
