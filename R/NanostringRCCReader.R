
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
#'@export`
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