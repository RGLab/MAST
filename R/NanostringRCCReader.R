
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
##' @rdname readSection
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
  key<-rename(read.csv(file,header=FALSE, as.is=TRUE),c("V1"="Name","V2"="GeneID"))
  return(lapply(rcc,function(x){x$code_summary<-merge(x$code_summary,key);x}))
}
 

setMethod('initialize', signature='NanoStringAssay', function(.Object, ...){
  .Object <- callNextMethod()
  if(max(exprs(.Object))>100) warning('log Counts > 100 found; are you sure the data was log2 + 1 transformed?')
  layername(.Object) <- 'lCount'
  .Object
})

##' Constructor for a NanoStringAssay
##'
##' Constructs a NanoStringAssay object. Differs little from the FluidigmAssay constructor. Accepts a function for post-processing.
##' mapping argument has been removed as well for simplicity.
##' @title NanoString Assay Constructor
##' @param rccfiles A list of rcc files with full paths
##' @param keyfile A full path to a keyfile
##' @param idvars See \code{\link{SingleCellAssay}}
##' @param primerid See \code{\link{SingleCellAssay}}
##' @param measurement The measurement column for raw data. 
##' @param ncells A \code{character} specifying the column which gives the number of cells per well
##' @param id An identifier for the resulting object. Should be a meaningful name
##' @param cellvars See \code{\link{SingleCellAssay}}
##' @param featurevars See \code{\link{SingleCellAssay}}
##' @param phenovars See \code{\link{SingleCellAssay}}
##' @param post.process.function function applied to \code{data.frame} of all rcc files, before the NanostringAssay object is constructed.
##' @param ... Additional parameters passed to \code{SingleCellAssay} constructor
##' @return A NanoStringAssay object
##' @author Andrew McDavid and Greg Finak
##' @export NanoStringAssay
NanoStringAssay<-function(rccfiles=NULL,keyfile=NULL,idvars,primerid,measurement='Count', ncells=NULL,id=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, post.process.function=NULL, ...){
  
  rcclist<-readNanoStringLanes(rccfiles);
  if(!is.null(keyfile)) rcclist<-mergeWithKeyFile(rcc=rcclist,keyfile)
  
  ## lengths <- sapply(rcclist, function(rccout){max(nrow(rccout$code_summary), nrow(rccout$lane_attributes), nrow(rccout$sample_attributes))})
  ## widths <- sapply(rcclist
  ## if(length(rcclist)>1) skeleton <- cbind(
                                               
                    
  suppressWarnings(
    dataframe<-rbindlist(lapply(1:length(rcclist),function(i){
    cbind(rcclist[[i]]$code_summary,rcclist[[i]]$lane_attributes,rcclist[[i]]$sample_attributes, stringsAsFactors=FALSE)
  }))
    )
  dataframe <- as.data.frame(dataframe)
 
  if(!is.null(post.process.function)){
    dataframe<-post.process.function(dataframe) 
  }
  dataframe$lCount <- log2(dataframe[,measurement]+1)
  measurement <- 'lCount'
  #reame the mapping 
  this.frame <- as.list(environment())
  #remove unneeded variables
  this.frame$post.process.function<-NULL
  this.frame$rcclist<-NULL
  this.frame$rccfiles<-NULL
  this.frame$keyfile<-NULL
  this.frame$keep.names <- FALSE
  this.frame$sort <- FALSE
  
  ## Factor out code that builds the mapping
  this.frame$cellvars <- c(ncells, cellvars)
  #sc <- do.call(SingleCellAssay, this.frame)
  #sc<-as(sc,'NanoStringAssay')
  cmap <- new('Mapping')
  cmap['ncells'] <- ncells
  sc <- new('NanoStringAssay', dataframe=dataframe, idvars=idvars, primerid=primerid, measurement=measurement, id=id, cellvars=cellvars, featurevars=featurevars, phenovars=phenovars, cmap=cmap)
  return(sc)
}
