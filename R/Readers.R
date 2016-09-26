##Function to read fluidigm data
## files = character vector of files (csv or xls)
## metadata: csv, maps to files
## header.size: currently ignored
## skip: currently ignored
## cycle.threshold: max number of cycles?
## metadataColClasses: used to infer datatypes in metadata (optional, then we guess datatypes)
## meta.key: ??
## idvars: vector of column names that uniquely identify a sample (conditional on splitby)
## splitby: values in this column are used to split input (per file??) into distinct SCA objects ala split
## unique.well.id: column name which identifies wells.  Must be unique within ???
## raw: assume files follow 'raw' schema
## assay:

##Function to read fluidigm data
##' Reads a fluidigm raw data file (or set of files)
##'
##' This function reads a raw Fluidigm data file or set of files and constructs a SingleCellAssay (or FluigidmAssay) object.
##' @title read.fluidigm
##' @name read.fluidigm
##' @param files A \code{character} vector of files to read.
##' @param metadata A \code{character} path and filename of a CSV file containing additional metadata about the samples
##' @param header.size A \code{numeric} indicating the number of lines in the header (default 2)
##' @param skip \code{numeric} how many lines to skip before reading (default 8)
##' @param cycle.threshold The maximum number of PCR cycles performed (default 40) \code{numeric}
##' @param metadataColClasses Optional \code{character} vector giving the column classes of the metadata file. See \link{read.table}.
##' @param meta.key Optional \code{character} vector that identifies the key column between the metadata and the fluidigm data
##' @param idvars Optional \code{character} vector that defines the set of columns uniquely identifying a well (unique cell, gene, and condition).
##' @param splitby Optional \code{character} that defines the column / variable used to split the resulting data into a list of SingleCellAssay, such that unique levels of \code{splitby} each fall into their own SingleCellAssay. Ususally the experimental unit subjected to different treatments.
##' @param unique.well.id The column that uniquely identifies a sample well in the data. Default is "Chamber.ID".
##' @param raw \code{logical} flag indicating this is raw data coming off the instrument. Thus we make some assumptions about the column names that are present. 
##' @param assay \code{character} name of a column that uniquely identifies an Assay (i.e. gene). Default is NULL
##' @param geneid \code{character} names of the column that identifies a gene. Default is "Assay.Name"
##' @param sample \code{character} name of a column that uniquely identifies a sample
##' @param well \code{character} name of a column that uniquely identifies a well. Default "Well".
##' @param measurement \code{character} name of the column that holds the measurement. Default "X40.Ct".
##' @param measurement.processed \code{character} one of "Ct","40-Ct", or "et". If not "Ct", the measurement will be transformed.
##' @param ncells The column with the number of cells in this well.
##' @return list of \code{SingleCellAssay} holding the data.
##' @export read.fluidigm
##' @author Greg Finak
##' @importFrom utils read.csv
read.fluidigm<-function(files=NULL,metadata=NULL,header.size=2,skip=8,cycle.threshold=40,metadataColClasses=NULL,meta.key=NULL,idvars=NULL,splitby=NULL,unique.well.id="Chamber.ID",raw=TRUE,assay=NULL,geneid="Assay.Name",sample=NULL,well="Well",measurement="X40.Ct",measurement.processed="Ct",ncells="SampleRConc"){
    measurement.processed<-match.arg(measurement.processed,c("Ct","40-Ct","et"))
    if(raw){
        geneid<-"Gene"
        measurement<-"et"
        unique.well.id="Chamber.ID"
    }
    if(!is.null(metadata)){
        ##case no metadata
        ##Read the metadata file
        message("\nReading Metadata")
        metadata<-read.csv(metadata,colClasses=metadataColClasses,header=TRUE)
    }else if(raw){
        ##case no metadat and raw data
        stop("Must provide metadata file when loading raw data.")
    }
    
    
    if(!any(colnames(metadata)%in%meta.key)&!is.null(meta.key)){
        stop("Column \"", meta.key, "\" not found in metadata file")
    }
    if(is.null(idvars)){
        stop("Must provide `idvars` specifying columns that identify unique rows conditional on `splitby`")
    }
    if(is.null(files)){
        stop("Invalid file name")
    }
    if(is.null(unique.well.id)){
        stop("unique.well.id cannot be NULL")
    }
    nfiles<-length(files)
    
    ##Read raw data files if raw is TRUE
    if(raw){
        message("\nReading ",nfiles," files.")
        if(any(grepl("parallel",loadedNamespaces()))){
            file<-as.list(files)
            datas<-parallel::mclapply(files,read.fluidigm.xls)
        }else{
            datas<-lapply(files,read.fluidigm.xls)
        }
    }else{
        ##read processed data.
        ##check if they're csv files
        if(!all(grepl("csv$",files))){
            stop("Processed data must be in .csv files")
        }
        files<-as.list(files)
        if(any(grepl("parallel",loadedNamespaces()))){
            datas<-parallel::mclapply(files,read.csv)
        }else{
            datas<-lapply(files,read.csv)                                               
        }
        ##Set the filename attribute for matching to metadata if necessary
        names(datas)<-files
        lapply(names(datas),function(x){attr(datas[[x]],"filename")<-x})        
    }                                   
    ##Match the filename using meta.key user specificed column.. if metadata is present and it's raw...
    ##case metadata and raw or processed data - match
    if(!is.null(metadata)){
        if(!is.null(meta.key)){
            form<-eval(substitute(~meta.key))
            form[[2]]<-as.name(form[[2]])
            j<-na.omit(match(as.matrix(model.frame(form,metadata)),basename(unlist(lapply(datas,function(x)attributes(x)$filename),use.names=FALSE))))
        }
    }
    
    
    if(raw){    
        if(!all(unlist(lapply(datas,function(x)unique.well.id%in%colnames(x)),use.names=FALSE))){
            stop("unique.well.id = ",unique.well.id, " is not valid for some of the data files")
        }
        for(i in 1:length(datas)){
            datas[[i]]<-within(datas[[i]],{
                sample<-factor(do.call(c,lapply(strsplit(as.character(get(unique.well.id)),"-"),function(x)x[1])))
                assay<-factor(do.call(c,lapply(strsplit(as.character(get(unique.well.id)),"-"),function(x)x[2])))
            }
            )
        }
    } 
    else if(is.null(well)&!is.null(geneid)&!is.null(sample)&!is.null(assay)) {
        ##TODO fill this in if needed
    } 
    else if(!is.null(well)&!is.null(geneid)) {
        ##TODO fill this in if needed.
    } 
    else {
        stop("You must provide the column mapping for the well, or the sample and assay, as well as the gene when parsing processed data")
    }
    
    if(!raw) {
        ##construct a primerid and et columns for processed data.
        for(i in seq_along(datas)){
            datas[[i]]<-within(datas[[i]],{
                primer.id<-get(geneid)
                switch(measurement.processed,"Ct"= {et<-40-get(measurement); et[et<0]<-0}
                      ,
                       "40-Ct"= {et<-get(measurement);et[et<0]<-0}
                      ,
                       "et"= {et<-get(measurement)}
                       )
            })
        }
    }
    
    
    
    ##construct unique primer IDs
    if(raw){
        for(i in seq_along(datas)){
            pid<-with(datas[[i]],{
                intervals<-apply(apply(apply(table(assay,Gene),2,function(x)x>0),2,cumsum),2,function(x)findInterval(c(0:10),x))+1
                nl<-length(levels(factor(assay)))
                gene.assay.map<-sapply(levels(Gene),function(g){
                    rownames(table(assay,Gene))[intervals[intervals[,g]<nl,g]]
                })
                probe.id<-as.character(Gene)
                for(g in names(gene.assay.map)){
                    if(length(gene.assay.map[[g]])==0){
                        ##gene is not on the array
                    }
                    else
                    {
                        for(j in seq_along(gene.assay.map[[g]])){
                            pname<-paste(probe.id[assay%in%gene.assay.map[[g]][[j]]],j,sep=".")
                            probe.id[assay%in%gene.assay.map[[g]][[j]]]<-pname
                        }       
                    }
                }
                return(probe.id)
            })
            datas[[i]]$primer.id<-factor(pid)
        }
        
        ##construct an et column
        for(i in seq_along(datas)){
            datas[[i]]<-within(datas[[i]],{
                et<-`Ct Value`
                et[et==999]<-cycle.threshold
                et<-40-et
            }
            )
        }
        for(i in seq_along(j)){
            a<-attributes(datas[[j[[i]]]])
            d<-data.frame(datas[[j[[i]]]],metadata[i,])
            a$names<-attr(d,"names")
            datas[[j[[i]]]]<-d
            attributes(datas[[j[[i]]]])<-a
        }
    }
    bound<-do.call(rbind,datas)
    ##construct a single cell assay
    ##need an expression matrix, a cell-level annotated data frame, a feature-level annotated data frame, and a phenotype-level annotated data frame
    ##expression matrix should be cells x genes.
    ##splitby, idvars
    ##constructor for SCA set
    ##constructor for single cell assay or other type of assay

    joint <- FromFlatDF(dataframe=bound, idvars=idvars,primerid="primer.id",measurement="et",featurevars=geneid,ncells=ncells,class="FluidigmAssay")
    if(is.null(splitby)) return(joint)
    set<-split(joint, splitby)
    return(set)
}



if(getRversion() >= "2.15.1") globalVariables(c('Ct Call', 'Ct Value', 'Gene'))


##Function to read fluidigm from xls file
read.fluidigm.xls<-function(x,header.size=2,skip=8){
    colmap<-c(`Chamber ID`="Chamber.ID",`Sample Name`="SampleName",`Sample Type`="SampleType",`Sample rConc`="SampleRConc",`FAM-MGB Name`="Gene",`FAM-MGB Type`="AssayType",`Ct Quality`="CtQuality",`Ct Threshold`="CtThreshold")
    if(!class(x)=="character"){
        message("Argument must be of type character")
    }
    if(!file.exists(x)){
        stop("File not found ",x)
    }
    if(grepl("xls$",basename(x))){
        ##Read as xls file
        tmp<-gdata::xls2csv(x)
    }else if(grepl("csv$",basename(x))){
        ##read as csv
        tmp<-file(x,open="r")
    }
    header<-readLines(tmp,8)
    quality.threshold<-strsplit(gsub("\\\"","",header[5]),",")[[1]][2]
    baseline.correction<-strsplit(gsub("\\\"","",header[6]),",")[[1]][2]
    
    filename<-basename(x);
    
    cnames<-readLines(tmp,10)[9:10]
    cnames<-gsub("\\\"","",do.call(paste,lapply(cnames,function(x)strsplit(x,",")[[1]])))
    data<-read.csv(tmp,skip=10,header=FALSE)
    colnames(data)<-cnames
    data<-rename(data,colmap)
    data<-subset(data,!(`Ct Call`%in%"Fail"&`Ct Value`<999))
    attr(data,"quality.threshold")<-quality.threshold
    attr(data,"baseline.correction")<-baseline.correction
    attr(data,"filename")<-filename             
    return(data)
}
