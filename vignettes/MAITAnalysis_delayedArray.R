## ----libraries,echo=FALSE, error=FALSE,echo=TRUE,results='hide',warning=FALSE----
suppressPackageStartupMessages({
    library(ggplot2)
    library(GGally)
    library(GSEABase)
    library(limma)
    library(reshape2)
    library(data.table)
    library(knitr)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(stringr)
    library(NMF)
    library(rsvd)
    library(RColorBrewer)
    library(MAST)
  library(HDF5Array)
  library(pbapply)
  library(mbenchmark)
})
#options(mc.cores = detectCores() - 1) #if you have multiple cores to spin
options(mc.cores = 1)
# knitr::opts_chunk$set(message = FALSE,error = FALSE,warning = FALSE,cache = FALSE,fig.width=8,fig.height=6, auto.dep=TRUE)

## ----MASToptions---------------------------------------------------------
freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

## ----data,results='hide'-------------------------------------------------
data(maits, package='MAST')
head(maits$cdat)
head(maits$fdat)
mat <- t(maits$expressionmat)
dim(mat)

#' Convert bigmemory-version of mat
bm <- as.big.matrix(mat, backingpath = tempdir(), backingfile = "bm.bin", descriptorfile = "bm.desc")
file.name(bm)
bm <- DelayedArray(BMArraySeed(bm))
rownames(bm) <- rownames(mat)
colnames(bm) <- colnames(mat)
all.equal(mat, as.matrix(bm))

## ----createSca-----------------------------------------------------------
scaRaw <- FromMatrix(mat, maits$cdat, maits$fdat)
scaRaw.bm <- FromMatrix(bm, maits$cdat, maits$fdat)

# scaRaw.bm <- scaRaw
# assay(scaRaw.bm) <- bm # hack

# all.equal(assay(scaRaw), mat)
# ----heatmap doesn't work basically `NMF` package isn't aware of `DelayedArray::rowMeans`, and question is if there is anything I can do about it without asking `NMF` author to import `BiocGenerics::rowMeans`?
# aheatmap(assay(scaRaw[1:1000,]), labRow='', annCol=as.data.frame(colData(scaRaw)[,c('condition', 'ourfilter')]), distfun='spearman')
# aheatmap(assay(scaRaw.bm[1:1000,]), labRow='', annCol=as.data.frame(colData(scaRaw.bm)[,c('condition', 'ourfilter')]), distfun='spearman')

## ----pca works out-of-box
set.seed(123)
plotPCA <- function(sca_obj){
    projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
    colnames(projection)=c("PC1","PC2","PC3","PC4")
    pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
    print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'libSize', 'PercentToHuman', 'nGeneOn', 'exonRate'),
            mapping=aes(color=condition), upper=list(continuous='blank')))
    invisible(pca)
}

# plotPCA(scaRaw)
# plotPCA(scaRaw.bm)

filterCrit <- with(colData(scaRaw), pastFastqc=="PASS"& exonRate >0.3 & PercentToHuman>0.6 & nGeneOn> 4000)

#' filtering
sca <- subset(scaRaw,filterCrit)
sca.bm <- subset(scaRaw.bm,filterCrit)

eid <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys = mcols(sca)$entrez,keytype ="GENEID",columns = c("GENEID","TXNAME"))
ueid <- unique(na.omit(eid)$GENEID)

sca <- sca[mcols(sca)$entrez %in% ueid,]
sca.bm <- sca.bm[mcols(sca.bm)$entrez %in% ueid,]

## Remove invariant genes
set.seed(123)
sca <- sca[sample(which(freq(sca)>0), 6000),]
set.seed(123)
sca.bm <- sca.bm[sample(which(freq(sca.bm)>0), 6000),]

all.equal(assay(sca), as.matrix(assay(sca.bm)))

## ---- fig.width=4, fig.height=4------------------------------------------
cdr2 <-colSums(assay(sca)>0)
cdr2.bm <-colSums(assay(sca.bm)>0)
all.equal(cdr2, cdr2.bm)
# qplot(x=cdr2, y=colData(sca)$nGeneOn) + xlab('New CDR') + ylab('Old CDR')

## ------------------------------------------------------------------------
colData(sca)$cngeneson <- scale(cdr2)
colData(sca.bm)$cngeneson <- scale(cdr2)
## ----pcaFilter, dpi = 36-------------------------------------------------
# plotPCA(sca.bm)

## ----distribution, fig.width=6, fig.height=6-----------------------------
set.seed(123)
scaSample <- sca[sample(which(freq(sca)>.1), 20),]
set.seed(123)
scaSample.bm <- sca.bm[sample(which(freq(sca.bm)>.1), 20),]
all.equal(assay(scaSample), as.matrix(assay(scaSample.bm)))

flat <- as(scaSample, 'data.table')
flat.bm <- as(scaSample.bm, 'data.table')
all.equal(flat, flat.bm)
# ggplot(flat, aes(x=value))+geom_density() +facet_wrap(~symbolid, scale='free_y')

## ----threshold, results='hide', fig.width=6------------------------------
setRealizationBackend("BMArray")

thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
thres.bm <- thresholdSCRNACountMatrix(assay(sca.bm), nbins = 20, min_per_bin = 30)
object.size(thres)
object.size(thres.bm)

all.equal(thres$counts_threshold, as.matrix(thres.bm$counts_threshold))
all.equal(thres$original_data, as.matrix(thres.bm$original_data))
all.equal(thres[-(1:2)], thres.bm[-(1:2)])


a <- array(0, dim=c(3, 3))
b <- array(1, dim = c(2,2))
a
b
a[1:2, 1:2] <- b
a
A <- DelayedArray(a)
idx <- DelayedArray(array(c(T,T,F,T,T,F,F,F,F), dim  = c(3,3)))
A[idx] <- DelayedArray(b)
tt <- tempfile()
h5write(a, tt, "data")
A <- HDF5Array(tt, "data")
B
par(mfrow=c(5,4))
plot(thres)

## ----assignThresh--------------------------------------------------------
assays(sca) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
expressed_genes <- freq(sca) > freq_expressed
sca <- sca[expressed_genes,]

## ----zlm-----------------------------------------------------------------
cond<-factor(colData(sca)$condition)
cond<-relevel(cond,"Unstim")
colData(sca)$condition<-cond
zlmCond <- zlm(~condition + cngeneson, sca)
# The following are equivalent
## lrt <- lrTest(zlm, "condition")
## lrt <- lrTest(zlm, CoefficientHypothesis('conditionStim'))

# This would test if 2*cngeneson=conditionStim
#  This is sheer nonsense biologically and statistically, but gives an example of the flexibility.
## lrt <- lrTest(zlm, Hypothesis('2*cngeneson-conditionStim'))

## ----zlmSummary----------------------------------------------------------
#only test the condition coefficient.
summaryCond <- summary(zlmCond, doLRT='conditionStim') 
#print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')

## ------------------------------------------------------------------------
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='conditionStim' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='conditionStim' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)

## ---- fig.width=8,fig.height=8, dpi = 36---------------------------------
entrez_to_plot <- fcHurdleSig[1:50,primerid]
symbols_to_plot <- fcHurdleSig[1:50,symbolid]
flat_dat <- as(sca[entrez_to_plot,], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=condition, y=thresh,color=condition)) + geom_jitter()+facet_wrap(~symbolid, scale='free_y')+ggtitle("DE Genes in Activated MAIT Cells")
ggbase+geom_violin() 

## ---- dpi = 36-----------------------------------------------------------
flat_dat[,lmPred:=lm(thresh~cngeneson + condition)$fitted, key=symbolid]
ggbase +aes(x=cngeneson) + geom_line(aes(y=lmPred), lty=1) + xlab('Standardized Cellular Detection Rate')


## ---- dpi = 36-----------------------------------------------------------
## This is all rather kludgy at the moment
MM <- model.matrix(~condition,unique(colData(sca)[,c("condition"),drop=FALSE]))
rownames(MM) <- str_extract(rownames(MM), 'Stim|Unstim')
predicted <- predict(zlmCond,modelmatrix=MM)

## Avert your eyes...
predicted[, primerid:=as.character(primerid)]
predicted_sig <- merge(mcols(sca), predicted[primerid%in%entrez_to_plot], by='primerid')
predicted_sig <- as.data.table(predicted_sig)

## plot with inverse logit transformed x-axis
ggplot(predicted_sig)+aes(x=invlogit(etaD),y=muC,xse=seD,yse=seC,col=sample)+
    facet_wrap(~symbolid,scales="free_y")+theme_linedraw()+
    geom_point(size=0.5)+scale_x_continuous("Proportion expression")+
    scale_y_continuous("Estimated Mean")+
    stat_ell(aes(x=etaD,y=muC),level=0.95, invert='x')


## ---- fig.height=8-------------------------------------------------------
mat_to_plot <- assay(sca[entrez_to_plot,])
rownames(mat_to_plot) <- symbols_to_plot
aheatmap(mat_to_plot,annCol=colData(sca)[,"condition"],main="DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)))


## ------------------------------------------------------------------------
table(colData(sca)$beta, exclude=NULL)
#Note that we currently throw an uninformative error if a covariate is `NA`
scaHasBeta <- subset(sca, !is.na(beta))

## ----residuals-----------------------------------------------------------
scaDE <- sca[entrez_to_plot,]
zlmResidDE <- zlm(~condition + cngeneson, scaDE, hook=deviance_residuals_hook)
residDE <- zlmResidDE@hookOut
residDEMatrix <- do.call(rbind, residDE)

## ----addResiduals, dpi = 36----------------------------------------------
assays(scaDE) <- c(assays(scaDE), list(resid=residDEMatrix))
scaResidFlat <- as(scaDE, 'data.table')
scaResidFlat[1:4,]
ggplot(scaResidFlat, aes(x=ngeneson, y=resid))+geom_point(aes(col=condition))+geom_smooth()+facet_wrap(~symbolid)


## ----boots, eval=TRUE----------------------------------------------------
# bootstrap, resampling cells
# R should be set to >50 if you were doing this for real.
boots <- bootVcov1(zlmCond, R=4)

## ----setupBTM, dependson="data;zlm;boots",results='hide'-----------------
module <- "BTM"
min_gene_in_module <- 5
packageExt <- system.file("extdata", package='MAST')
module_file <- list.files(packageExt, pattern = module, full.names = TRUE)
gene_set <- getGmt(module_file)
gene_ids <- geneIds(gene_set)
gene_ids <- gene_ids[!names(gene_ids)%like%"TBA"&!names(gene_ids)%like%"B cell"]
sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$symbolid)
# Only keep modules with at least min_gene_in_module
sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]


## ----gsea----------------------------------------------------------------
gsea <- gseaAfterBoot(zlmCond, boots, sets_indices, CoefficientHypothesis("conditionStim")) 
z_stat_comb <- summary(gsea, testType='normal')

## ----gseaView------------------------------------------------------------
sigModules <- z_stat_comb[combined_adj<.01]
gseaTable <- melt(sigModules[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')
ggplot(gseaTable, aes(y=set, x=variable, fill=value))+geom_raster() + scale_fill_distiller(palette="PiYG")

