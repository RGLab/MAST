#' Makes a nice BiPlot
#'
#' Creates a custom BiPlot for visualizing the results of PCA
#' @param pc output of \code{prcomp}
#' @param colorfactor a \code{factor} the same length as \code{nrow(pc$x)} to color the points
#' @param scaling \code{integer} to scale the vectors showing loadings
#' @param nudge \code{numeric} to offset labels for loadings
#' @param N number of variables with longest \code{dim[1]} or \code{dim[2]} projections to display
#' @param dims \code{numeric} vector of length 2 indicating which PCs to plot
#' @param ... passed to plot
#' @return printed plot
myBiplot <- function(pc,colorfactor,scaling=100,nudge=1.2,N=10,dims=1:2,...){
    o <- unique(c(order(abs(pc$rotation[,dims[1]]),decreasing=TRUE)[1:N],order(abs(pc$rotation[,dims[2]]),decreasing=TRUE)[1:N]))
    colors <- colorfactor
    textNudge <- nudge
    plot(pc$x[,dims[1]],pc$x[,dims[2]],col=colors,asp=1,las=1,type="n",xlab=paste0("PC",dims[1]),ylab=paste0("PC",dims[2]),...)
    arrows(0,0,pc$rotation[o,dims[1]]*scaling,scaling*pc$rotation[o,dims[2]],col="blue",length=0.1)
    text(pc$x[,dims[1]], pc$x[,dims[2]], rownames(pc$x), col=as.integer(colors), cex=0.7)
    text(pc$rotation[o,dims[1]]*scaling*textNudge, textNudge*scaling*pc$rotation[o,dims[2]], rownames(pc$rotation[o,]), col="purple", cex=0.7)
}
