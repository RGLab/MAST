library(SingleCellAssay)

#develop a reader for nanostring data

datafiles<-list.files(pattern="RCC",path="/Users/gfinak/Documents/Projects/Nanostring/H9Cells/",recursive=TRUE,full=TRUE)
keyfiles<-"/Users/gfinak/Documents/Projects/Nanostring/h9key.csv"







process<-function(full=NULL){
#Date
full$Date<-as.Date(as.character(full$Date),format="%Y%m%d")

#Plate to Cell Cycle mapping
full<-merge(full,data.frame(Plate=factor(c(1,2,3,4,5,6)),CellCycle=c("G0/G1","G0/G1","S","S","G2/M","G2/M")))

#Number of cells is constant (1 cell per lane)
full$ncells<-1

#parse the comments in the sample info to get the cell type, plate, and row ID
full<-cbind(full,data.frame(CellType=substring(as.character(full$Comments),1,2),
                            Plate=substring(as.character(full$Comments),3,3),
                            RowID=substring(as.character(full$Comments),4,4)))
full
}
F<-NanoStringAssay(rccfiles=datafiles[1:10],keyfile=keyfiles,idvars=c("Plate","CartridgeID","ID"),geneid="GeneID",measurement="Count",phenoVars=c("CellType","CellCycle"),ncells="ncells",primerid="GeneID",cellvars=c("ncells","CellType"),id="nanostring",post.process.function=process)

