
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
    con<-file(x[i],open="r")#open read only
    lanes[[i]]$header<-.readHeader(con)
    lanes[[i]]$sample_attributes<-.readSampleAttributes(con)
    lanes[[i]]$lane_attributes<-.readLaneAttributes(con)
    lanes[[i]]$code_summary<-.readCodeSummary(con)
    close(con)
  }
  return(lanes)
}

#'Read the RCC header
#'
#'@param x is a connection
#'@return a data frame
.readHeader<-function(x){
  STARTHEADER<-"<Header>"
  ENDHEADER<-"</Header>"
  curline<-readLines(x,n=1)
  if(curline!=STARTHEADER){
    desc<-summary(con)$description
    close(con)
    stop("Oops.. failed reading RCC file header for ",desc)
  }
  repeat{
    curline<-c(curline,readLines(x,n=1))
    if(curline[length(curline)]==ENDHEADER){
      break
    } 
  }
  header<-read.csv(textConnection(paste(curline[-c(1,length(curline))],collapse="\n")),header=FALSE)
  return(cast(melt(header,id="V1"),variable~V1)[,-1L])
}


#'Read the RCC Sample Attributes
#'
#'@param x is a connection
#'@return a data frame
.readSampleAttributes<-function(x){
  START<-"<Sample_Attributes>"
  END<-"</Sample_Attributes>"
  repeat{
    curline<-readLines(x,n=1)
    if(curline!=""){
      break
    }
  }
  if(curline!=START){
    desc<-summary(con)$description
    close(con)
    stop("Oops.. failed reading RCC file Sample Attributes for ",desc)
  }
  repeat{
    curline<-c(curline,readLines(x,n=1))
    if(curline[length(curline)]==END){
      break
    } 
  }
  sample_attributes<-read.csv(textConnection(paste(curline[-c(1,length(curline))],collapse="\n")),header=FALSE)
  return(cast(melt(sample_attributes,id="V1"),variable~V1)[,-1L])
}


#'Read the RCC Lane Attributes
#'
#'@param x is a connection
#'@return a data frame
.readLaneAttributes<-function(x){
  START<-"<Lane_Attributes>"
  END<-"</Lane_Attributes>"
  repeat{
    curline<-readLines(x,n=1)
    if(curline!=""){
      break
    }
  }
  if(curline!=START){
    desc<-summary(con)$description
    close(con)
    stop("Oops.. failed reading RCC file Sample Attributes for ",desc)
  }
  repeat{
    curline<-c(curline,readLines(x,n=1))
    if(curline[length(curline)]==END){
      break
    } 
  }
  lane_attributes<-read.csv(textConnection(paste(curline[-c(1,length(curline))],collapse="\n")),header=FALSE)
  return(cast(melt(lane_attributes,id="V1"),variable~V1)[,-1L])
}

#'Read the RCC Lane Attributes
#'
#'@param x is a connection
#'@return a data frame
.readCodeSummary<-function(x){
  START<-"<Code_Summary>"
  END<-"</Code_Summary>"
  repeat{
    curline<-readLines(x,n=1)
    if(curline!=""){
      break
    }
  }
  if(curline!=START){
    desc<-summary(con)$description
    close(con)
    stop("Oops.. failed reading RCC file Sample Attributes for ",desc)
  }
  repeat{
    curline<-c(curline,readLines(x,n=1))
    if(curline[length(curline)]==END){
      break
    } 
  }
  code_summary<-read.csv(textConnection(paste(curline[-c(1,length(curline))],collapse="\n")),header=TRUE)
  return(code_summary)
}

#'Merge NanoString lanes with a key file
#'
#'Reads in a key file and maps to nanostring lanes
#'@param rcc is a list returned from readNanoStringLanes
#'@param file is a character vector name of the key file
#'@return mapped rcc list
#'@export
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
  dataframe<-do.call(rbind,lapply(1:length(rcclist),function(i){
    cbind(rcclist[[i]]$code_summary,rcclist[[i]]$lane_attributes,rcclist[[i]]$sample_attributes)  
  }))
  if(!is.null(post.process.function)){
    dataframe<-post.process.function(dataframe) 
  }
  #reconstruct this for the next call
  #transform with log+1
  dataframe$lCount<-log(get(measurement,dataframe)+1)
  
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
