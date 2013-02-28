
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
setClass('NanoStringAssay', contains='FluidigmAssay', prototype=prototype(mapNames=SingleCellAssay:::FluidigmMapNames,mapping=new("Mapping",mapping=SingleCellAssay:::FluidigmMap)),validity=SingleCellAssayValidity)


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
##' @param geneid See \code{\link{SingleCellAssay}}
##' @param id An identifier for the resulting object. Should be a meaningful name
##' @param cellvars See \code{\link{SingleCellAssay}}
##' @param featurevars See \code{\link{SingleCellAssay}}
##' @param phenovars See \code{\link{SingleCellAssay}}
##' @param post.process.function function applied to \code{data.frame} of all rcc files, before the NanostringAssay object is constructed.
##' @param ... Additional parameters passed to \code{SingleCellAssay} constructor
##' @return A FluidigmAssay object
##' @author Andrew McDavid and Greg Finak
##' @export NanoStringAssay
NanoStringAssay<-function(rccfiles=NULL,keyfile=NULL,idvars,primerid,measurement, ncells=NULL, geneid=NULL,id=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, post.process.function=NULL,...){
  mapping<-new("Mapping",mapping=SingleCellAssay:::FluidigmMap)
  
  #transform the counts and rename the mapping for measurement to the new name
  #TODO fix the duplication of column names in a more reasonable way
  if(!is.null(ncells)){
    mapping<-addMapping(mapping,list(ncells=ncells))
  }
  
  rcclist<-readNanoStringLanes(rccfiles);
  rcclist<-mergeWithKeyFile(rcc=rcclist,keyfile)
  suppressWarnings(
    dataframe<-do.call(rbind,lapply(1:length(rcclist),function(i){
    cbind(rcclist[[i]]$code_summary,rcclist[[i]]$lane_attributes,rcclist[[i]]$sample_attributes, stringsAsFactors=FALSE)  
  }))
    )
  ## dataframe <- lapply(dataframe, function(col){
  ##   suppressWarnings(numTry <- as.numeric(col))
  ##   ## non-numeric character vectors return NA
  ##   ## if all entries that weren't originally NA are not NA, then we suppose everything was a numeric
  ##   if( all(!is.na(numTry[!is.na(col)])) ){
  ##     return(numTry)
  ##   }
  ##   return(col)
  ## })
  if(!is.null(post.process.function)){
    dataframe<-post.process.function(dataframe) 
  }
  #reconstruct this for the next call
  #transform with log+1
  dataframe$lCount<-log2(get(measurement,dataframe)+1)
  
  #reame the mapping 
  raw<-measurement
  measurement<-"lCount"
  this.frame <- as.list(environment())
  #remove unneeded variables
  this.frame$post.process.function<-NULL
  this.frame$rcclist<-NULL
  this.frame$rccfiles<-NULL
  this.frame$keyfile<-NULL
  
  ## Factor out code that builds the mapping
  this.frame$cellvars <- c(ncells, cellvars)
  sc <- do.call(SingleCellAssay, this.frame)
  sc@mapping<-addMapping(getMapping(sc),list(raw="Count"))
  #TODO need a cast for NanoStringAssay
  sc<-as(sc, 'FluidigmAssay')
  sc<-as(sc,'NanoStringAssay')
  return(sc)
}

setAs('FluidigmAssay', 'NanoStringAssay', function(from)  new("NanoStringAssay",env=from@env,mapping=addMapping(from@mapping,list(ncells=NULL)),id=from@id, wellKey=from@wellKey, featureData=from@featureData, phenoData=from@phenoData, cellData=from@cellData, description=from@description))
