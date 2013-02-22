library(SingleCellAssay)

#develop a reader for nanostring data

datafiles<-list.files(pattern="RCC",path="/Users/gfinak/Documents/Projects/Nanostring/H9Cells/",recursive=TRUE,full=TRUE)
keyfiles<-"/Users/gfinak/Documents/Projects/Nanostring/h9key.csv"


key<-read.csv(keyfiles)


foo<-readNanoStringLanes(datafiles[1:20])
