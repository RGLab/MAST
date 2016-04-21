##' Filter an RNASeqAssay
##'
##' Use the percentage of library per gene CDF type plots
##' to pick out outliers using the \code{1.5*IQR}
##' threshold
##' @param assay a \code{RNASeqAssay}
##' @param plot a \code{logical} indicating whether plots should be generated
##' @export
##' @return a filtered RNASeqAssay
outlierDetection <- function(assay,plot=TRUE){
  if(!is(assay,"RNASeqAssay")){
    stop(sprintf("assay must be an RNASeqAssay, but found %s",class(assay)[1]))
  }
  
  #compute the %of total library size for each gene
  f <- (exprs(assay)) / (rowSums(exprs(assay)))
  
  #filter out unexpressed genes
  gfilt <- apply(f,2,sum)!=0
  f <- f[,gfilt]
  
  #cumulative sum by library
  fc <- t(apply(f,1,cumsum))
  
  #median
  md <- apply(fc,1,median)
  
  #filter at 1.5*IQR
  filter <- (md>(median(md)+1.5*IQR(md))|md<(median(md)-1.5*IQR(md)))
  
  if(plot){
    
    #'plot the outliers
    par(mfrow=c(1,2))
    plot(y=md,x=rank(md),main="Outliers",xlab="rank",ylab="median frequency")
    points(y=md[filter],x=rank(md)[filter],col="red")
    
    #cdfs
    plot(fc[1,],type="l",main="Outliers",xlab="rank",ylab="cdf")
    apply(fc,1,lines)
    apply(fc[filter,],1,function(x)lines(x,col="red"))
  }
  #return the filtered object
  return(assay[!filter,gfilt])
}


#'Filter low-expressing genes
#'
#'Filter out genes that have less than some percent threshold expression across all libraries
#'
#'@param assay a \code{SingleCellAssay} object
#'@param threshold a \code{numeric} between 0, and 1, specifying the threshold frequency below which genes will be filtered out
#'@export
filterLowExpressedGenes<-function(assay,threshold=0.1){
  if(!is(assay,"SingleCellAssay")){
    stop(sprintf("assay must be a SingleCellAssay but found %s\n",class(assay)[1]))
  }
  e<-exprs(assay)
  gp<-colMeans(e>0)
  pfilt<-ecdf(gp)(threshold)
  cat(sprintf("Filtered out %2.2f percent of genes at %2.2f percent threshold\n",pfilt*100,threshold*100))
  return(assay[,gp>threshold])
}

