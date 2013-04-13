
#'RCC Nanostring Reader
#'
#'Reads RCC Nanostring Lane files
#'and constructs a SingleCellAssay object
#'
#'Open each lane file and read line by line constructing the different components of the file
#'
#'Metadata is present in the files. One lane per file.
#'@param x is a \code{character} vector of full path names to RCC files
#'@return a list
#' @importFrom stringr str_trim
#'@export
readNanoStringLanes<-function(x){
  nfiles<-length(x)
  lanes<-vector("list",nfiles)
  for(i in 1:length(x)){
    this.file <- x[i]
    con<-str_trim(readLines(this.file))
    lanes[[i]]$header<-.readSection(con, this.file, 'Header', header=FALSE, as.is=TRUE)
    lanes[[i]]$sample_attributes<-.readSection(con, this.file, 'Sample_Attributes', header=FALSE, as.is=TRUE)
    lanes[[i]]$lane_attributes<-.readSection(con, this.file, 'Lane_Attributes', header=FALSE, as.is=TRUE)
    lanes[[i]]$code_summary<-.readSection(con, this.file, 'Code_Summary', recast=FALSE, header=TRUE, as.is=TRUE)
  }
  return(lanes)
}



##' Read a section of a Nanostring RCC file
##'
##' The section is surrounded by <type> </type>
##' @param x character vector of lines from RCC file, newlines stripped
##' @param file name 
##' @param type the section of the file to be read in
##' @param recast should the section be transposed?
##' @param ... argument spassed to read.csv
##' @return data.frame of section
.readSection <- function(x, file, type, recast=TRUE, ...){
START <- sprintf('<%s>', type)
END <- sprintf('</%s>', type)
si <- match(START, x)
ei <- match(END, x)
if(is.na(si) || is.na(ei) || ei<si){
  stop('Malformed ', type, ' in ', file)
}
  section<-read.csv(textConnection(x[(si+1):(ei-1)]),...)
if(recast){
  wide <- data.frame(t(section), stringsAsFactors=FALSE)[-1,]
  names(wide) <- t(section)[1,]
} else{
  wide <- section
}
 wide
}

#'Merge NanoString lanes with a key file
#'
#'Reads in a key file and maps to nanostring lanes
#'@param rcc is a list returned from readNanoStringLanes
#'@param file is a character vector name of the key file
#'@return mapped rcc list
#'@export
#' @importFrom reshape rename
mergeWithKeyFile<-function(rcc,file){
  key<-rename(read.csv(file,header=FALSE),c("V1"="Name","V2"="GeneID"))
  return(lapply(rcc,function(x){x$code_summary<-merge(x$code_summary,key);x}))
}
 

#Could write a constructor that takes a post-processing function...
##' @exportClass NanoStringAssay
setClass('NanoStringAssay', contains='FluidigmAssay',validity=SingleCellAssayValidity)


##' Constructor for a NanoStringAssay
##'
##' Constructs a NanoStringAssay object. Differs little from the FluidigmAssay constructor. Accepts a function for post-processing.
##' mapping argument has been removed as well for simplicity.
##' @title NanoString Assay Constructor
##' @param rccfiles A list of rcc files with full paths
##' @param keyfile A full path to a keyfile
##' @param idvars See \code{\link{SingleCellAssay}}
##' @param primerid See \code{\link{SingleCellAssay}}
##' @param measurement The measurement column for raw data. Will be placed in a variable named \code{raw} in the environment storing the data. The transformed and thresholded data will be placed in data.
##' @param ncells A \code{character} specifying the column which gives the number of cells per well
##' @param id An identifier for the resulting object. Should be a meaningful name
##' @param cellvars See \code{\link{SingleCellAssay}}
##' @param featurevars See \code{\link{SingleCellAssay}}
##' @param phenovars See \code{\link{SingleCellAssay}}
##' @param post.process.function function applied to \code{data.frame} of all rcc files, before the NanostringAssay object is constructed.
##' @param ltrans Should the counts be log2 + 1 transformed?
##' @param ... Additional parameters passed to \code{SingleCellAssay} constructor
##' @return A FluidigmAssay object
##' @author Andrew McDavid and Greg Finak
##' @export NanoStringAssay
NanoStringAssay<-function(rccfiles=NULL,keyfile=NULL,idvars,primerid,measurement, ncells=NULL,id=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, post.process.function=NULL,ltrans=FALSE, ...){
  
  rcclist<-readNanoStringLanes(rccfiles);
  if(!is.null(keyfile)) rcclist<-mergeWithKeyFile(rcc=rcclist,keyfile)
  suppressWarnings(
    dataframe<-do.call(rbind,lapply(1:length(rcclist),function(i){
    cbind(rcclist[[i]]$code_summary,rcclist[[i]]$lane_attributes,rcclist[[i]]$sample_attributes, stringsAsFactors=FALSE)  
  }))
    )
 
  if(!is.null(post.process.function)){
    dataframe<-post.process.function(dataframe) 
  }
  #reconstruct this for the next call
  #transform with log+1
  if(ltrans){
  dataframe$lCount<-log2(get(measurement,dataframe)+1)
  measurement<-"lCount"
}
  
  #reame the mapping 
  this.frame <- as.list(environment())
  #remove unneeded variables
  this.frame$post.process.function<-NULL
  this.frame$rcclist<-NULL
  this.frame$rccfiles<-NULL
  this.frame$keyfile<-NULL
  
  ## Factor out code that builds the mapping
  this.frame$cellvars <- c(ncells, cellvars)
  sc <- do.call(FluidigmAssay, this.frame)
  sc<-as(sc,'NanoStringAssay')
  return(sc)
}

#' Estimate thresholds for positive expression
#'
#' Estimates per-gene x unit thresholds for positive expression and truncates values below this threshold
#' Uncertain values (in terms of posterior probability of membership) can be set to NA or rounded left or right
#' Thresholds are estimated using a Gaussian mixture model with prior supplied by population estimates.
#' @name threshold
#' @param nsa NanostringAssay object
#' @param groups groups to apply thresholding
#' @param thresholds data.frame of thresholds  of groups x genes user may specify
#' @param posteriorprob Min posterior probability of cluster membership for an observation to be truncated 
#' @param clip Should values with uncertain posterior probs be clipped 
#' @return modifies nsa in place
NULL
## setMethod('threshold', signature='NanostringAssay', function(nsa, groups, thresholds, posteriorprob, clip=c('left', 'right', 'NA'){
##   exprs.all <- log2(melt(nsa)[, getMapping(nsa, 'raw')]+1)
##   out <- Mclust(exprs.all, G=1:3, modelNames='V')
##   if(out$G != 2)
##     stop("Uhoh, didn't find two clusters in complete data set")
##   means <- out$param$mean
##   scales <- out$param$variance$scale

##   exprs.split <- split(exprs.all, melt(nsa)[, c(groups, getMapping(nsa, 'primerid'))])

##   lapp <- lapply(exprs.split, function(x){
##     out <- Mclust(x, G=1:2, modelNames='V', prior=priorControl(shrinkage=.03, mean=means, scale=scales))
##     print(out$G)
##     if(out$G==1){
##       print(out$param$mean)
##     }
##     out
##   })
  
## })
