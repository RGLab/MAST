## source('common-fixtures.R')


## geneid="Gene"
## primerid='Gene'
## measurement='et'
## idvars=c('Subject.ID', 'Chip.Number', 'Stim.Condition', 'Population', 'Well')
## ncells <- 'Number.of.Cells'
## phenovars=NULL
## cellvars='Experiment.Number'
## featurevars=NULL

## data(vbeta)
## test_that("vbeta can be loaded",{
##   expect_that(vbeta,is_a("data.frame"))
## })



context("Generating a complete and incomplete subset")
dat_complete <- vbeta
countComplete <- table(do.call(paste, dat_complete[,idvars]))
expect_that(all(countComplete==countComplete[1]), is_true())
dat_incomplete <- dat_complete[-seq(1, nrow(dat_complete), by=101),]
counts <- table(do.call(paste, dat_incomplete[,idvars]))
expect_that(all(counts == counts[1]), is_false())

blank <- dat_complete[1,]
## blankinst <- new('SingleCellAssay', dataframe=blank, idvars=idvars, primerid=geneid, measurement=measurement)
## blankinst.fd <- FluidigmAssay(dataframe=blank, idvars=idvars, primerid=geneid, measurement=measurement, ncells=ncells)
## test_that("Can create empty instance with pre-existing wellKey",{
##   expect_that(blankinst, is_a("SingleCellAssay"))
## })

## test_that("Can get wellKey",{
##   expect_that(getwellKey(blankinst)[[1]], equals(digest(paste(blankinst@env$data[,getMapping(blankinst@mapping)$idvars],collapse=" "))))
## })

vbeta$et <- ifelse(is.na(vbeta$Ct), 0, 40-vbeta$Ct)
fd <- FromFlatDF(vbeta, idvars=idvars, primerid=primerid, measurement=measurement, ncells=ncells, geneid=geneid, sort=TRUE)
test_that('could construct from flattened data.frame', {
    expect_that(fd, is_a('SingleCellAssay'))
})

test_that('Has dimnames', {
    expect_is(dimnames(fd)[[1]], 'character')
    expect_is(dimnames(fd)[[2]], 'character')
})

context('Test subsetting')
test_that('Subset columns by index, name, boolean', {
    asubset <- c(5, 1, 4, 10, 15)
    ss <- fd[,asubset]
    expect_equal(mcols(ss), mcols(fd))
    expect_equal(colData(ss), colData(fd)[asubset,])

    asubset <- c('Sub01 1 Stim(SEB) CD154+VbetaResponsive A07', 'Sub01 1 Stim(SEB) CD154+VbetaResponsive A01')
    ss <- fd[,asubset]
    expect_equal(mcols(ss), mcols(fd))
    expect_equal(colData(ss), colData(fd)[asubset,])
    
    asubset <- rep(c(TRUE, FALSE, TRUE, FALSE), length.out=dim(fd)[[2]])
    ss <- fd[,asubset]
    expect_equal(mcols(ss), mcols(fd))
    expect_equal(colData(ss), colData(fd)[asubset,])

})

test_that('Subset rows by index, name, boolean', {
    asubset <- c(5, 1, 4, 10, 15)
    ss <- fd[asubset,]
    expect_equal(mcols(ss, use.names=TRUE), mcols(fd, use.names=TRUE)[asubset,])
    expect_equal(colData(ss), colData(fd))

    asubset <- c('BAX', 'CCL2')
    ss <- fd[asubset,]
    expect_equal(mcols(ss, use.names=TRUE), mcols(fd, use.names=TRUE)[asubset,])
    expect_equal(colData(ss), colData(fd))
    
    asubset <- rep(c(TRUE, FALSE, TRUE, FALSE), length.out=dim(fd)[[1]])
    ss <- fd[asubset,]
    expect_equal(mcols(ss, use.names=TRUE), mcols(fd, use.names=TRUE)[asubset,])
    expect_equal(colData(ss), colData(fd))
})
    
    
test_that('Cell data and feature data are correctly assigned on construction', {
    vb.manip <- within(vbeta, {
        et[Stim.Condition=='Stim(SEB)'] <- 2000
        et[Stim.Condition!='Stim(SEB)' & Gene=='TGFB1'] <- 0
    })
    vb.manip <- vb.manip[sample(nrow(vb.manip)),]
    fd.manip <- FromFlatDF(vb.manip, idvars=c("Subject.ID", "Chip.Number", "Well"), primerid='Gene', measurement='et', ncells='Number.of.Cells', geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'), phenovars=c('Stim.Condition','Time'), id='vbeta all')
    expect_true(all(assay(subset(fd.manip, Stim.Condition=='Stim(SEB)'))==2000))
})

sc <- fd
test_that("Can load complete data", {
  tab <- table(melt.SingleCellAssay(sc)$wellKey)
  expect_that(tab, is_equivalent_to(countComplete))
})

test_that("Cellkey unique identifies a cell", {
  tab <- table(melt(sc)$wellKey, do.call(paste, melt(sc)[, idvars, with=FALSE]))
  expect_true(all(tab %in% c(0,75)))
  
})


context('test construction helper funcs')
suppressPackageStartupMessages(library(data.table))
  naframe <- data.table(var=rep(c(1, 2), each=3), na=c(NA, -9, NA, -9, NA, -9))
test_that("uniqueModNA doesn't include NA", {
    setkeyv(naframe, colnames(naframe))
  expect_equal(nrow(MAST:::uniqueModNA(naframe, include='var')), 2)
  expect_equal(nrow(MAST:::uniqueModNA(naframe[,-2, with=FALSE], include='var')), 2)
})

test_that('uniqueModNA works on multiple columns', {
    ## Now should return every row, since every row is unique
    naframe$extra <- 1:nrow(naframe)
    setkeyv(naframe, colnames(naframe))
    expect_equivalent(unique(naframe), MAST:::uniqueModNA(naframe, include='var'))
})

sci<- SingleCellAssay(dat_incomplete, idvars=idvars, primerid=geneid, measurement=measurement)
test_that("Completes incomplete data", {
  expect_equal(nrow(melt(sci)), nrow(dat_complete))

  incomplete <- rbind(melt(fd[1:20,1:20], value.name=measurement),
                      melt(fd[21:50, 11:30], value.name=measurement)) #equally sized primerid blocks
  fd.incomplete <- FromFlatDF(incomplete, idvars=idvars, primerid=primerid, measurement=measurement, ncells='ncells', geneid=geneid, keep.names=TRUE)
  expect_message(FromFlatDF(incomplete, idvars=idvars, primerid=primerid, measurement=measurement, ncells='ncells', geneid=geneid, keep.names=TRUE), 'incomplete')
  expect_equal(nrow(fd.incomplete), 50)
  expect_equal(ncol(fd.incomplete), 30)
  expect_true(any(is.na(assay(fd.incomplete))))
  expect_equal(assay(fd.incomplete[21:50, 11:30]),assay(fd[21:50, 11:30]))
  expect_equal(assay(fd.incomplete[1:20, 1:20]),assay(fd[1:20, 1:20]))
})

## No more mapping, hurray!
## context('testing cell and feature dictionaries')

## complete<-(MAST:::melt(sc))[,setdiff(colnames(MAST:::melt(sc)),"__wellKey")]
## scd <- SingleCellAssay(complete, mapping=getMapping(sc))

## test_that('Cell data has correct number of row/columns', {
##   expect_that(cData(scd), is_a('data.frame'))
##   expect_equivalent(nrow(scd), nrow(cData(scd)))
##   expect_equal(unique(cData(scd)), cData(scd))
##   expect_equivalent(ncol(cData(scd)), length(unique(c(phenovars, idvars, cellvars))))
## })
## scd <- new('SingleCellAssay', .Data=sc@.Data, cellData=cellData(sc), featureData=featureData(sc))
## test_that('Feature data has correct number of row/columns', {
##  expect_that(fData(scd), is_a('data.frame'))
##   expect_that(featureData(scd), is_a('AnnotatedDataFrame'))
##    expect_equivalent(ncol(scd), nrow(fData(scd)))
##   expect_equal(unique(fData(scd)), fData(scd))
##  #one extra column for the primerid
##   expect_equivalent(ncol(fData(scd)), length(unique(c(featurevars, geneid, primerid)))+1) 
## })

context("Testing methods")

## This makes more sense to me than to propagate new wells/features consisting entirely of NA
test_that('NAs throw error when subsetting', {
  expect_error(sc[, c(1:4, NA)])
  expect_error(sc[c(boolind, NA), ])
})

test_that("Throw error when indexing with factors", {
    expect_error(sc[, factor(c('B3GAT1'))])
    expect_error(sc[factor('A'),])
})


exprComplete <- assay(sc)
test_that("Exprs works", {
  measurement <- 'et'                #fix so melt renames column correctly
  expect_is(exprComplete, "matrix")
  expect_equal(nrow(exprComplete), nrow(sc))
  ind <- seq(1, nrow(dat_complete), by=1042)
  expect_equal(melt(sc)[ind,value], as.vector(exprComplete)[ind])
  geneandrow <- melt(sc)[1054,c(geneid, "wellKey"), with=FALSE]  
  thect <- melt(sc)[1054, value]
  expect_equivalent(exprComplete[ geneandrow[[1]], geneandrow[[2]]], thect)
})

test_that('Subset with TRUE is unchanged', {
  suball <- subset(sc, TRUE)
  expect_equal(suball, sc)
})

test_that('Subset with FALSE returns empty set', {
  subnone <- subset(sc, FALSE)
  expect_that(all.equal(sc, subnone), is_a('character'))
  expect_equal(ncol(subnone), 0)
})


test_that('Subset with names from SingleCellAssay works', {
    stim <- table(colData(sc)$Stim.Condition)[1]
  sub1 <- subset(sc, Stim.Condition == names(stim))
  expect_equivalent(ncol(sub1), stim)
})

test_that('Subset throws an intelligent error if thesubset cannot be evaluated', {
 expect_that(subset(sc, NOTPRESENT==fdsfjkl), throws_error('not found'))
})

context("SCASet works")
## test_that('Can construct', {
## aset <- SCASet(melt(sc), primerid=primerid, idvars=idvars, measurement='et', splitby='Subject.ID')
## })

test_that('Can split',{
        splat <- split(sc, cData(sc)$Subject.ID)
        expect_that(splat, is_a('list'))
        expect_equal(nrow(sc), nrow(splat[[1]]))
        expect_equal(ncol(sc), sum(sapply(splat, ncol)))
        splat.byfieldname <- split(sc, 'Subject.ID')
        expect_equal(splat.byfieldname, splat)
        splat <- split(sc, c('Subject.ID', 'Population'))
        expect_that(splat, is_a('list'))
        splat <- split(sc, list(factor(cData(sc)$Subject.ID), factor(cData(sc)$Population)))
        expect_that(splat, is_a('list'))
  
})

test_that('Replace works', {
    assay(sc, 1) <- matrix(-111, nrow=nrow(sc), ncol=ncol(sc))
    expect_true(all(assay(sc)==(-111)))
})

context("Backwards compatibility")

test_that('Exprs', {
    expect_equal(assay(sc), t(exprs(sc)))
    exprs(sc)[1, 10] <- -5
    expect_equal(assay(sc)[10, 1], -5)
})


context('Combine works')
doubleid <- data.frame(id1=c(1, 1, 2), id2=c(1, 2, 3), et=rep(3, 3), f1=rep('A', 3))
smallsc <- FromFlatDF(doubleid, idvars=c('id1', 'id2'), primerid='f1', measurement='et', id='1')

test_that('combine works', {
    spl <- split(smallsc, 'id1')
    c1 <- combine(spl[[1]], spl[[2]])
    expect_that(c1, is_a('SingleCellAssay'))
    expect_equal(ncol(c1), 3)
    c2 <- combine(spl[[1]], spl[[1]], spl[[2]])
    expect_that(c2, is_a('SingleCellAssay'))
})

test_that('combine throws error for non-conforming',{
  expect_error(combine(fd, spl[[1]]))
})

context('Test replace methods')


test_that('Can replace cData', {
    cDat <- colData(fd)
    cDat$foo <- "bar"
    colData(fd) <- cDat
    expect_true('foo' %in% names(colData(fd)))

    empty <- data.frame()
    expect_error(colData(fd) <- empty, 'DataFrame')

    scramble <- cDat[sample(nrow(cDat)),]
    expect_error(colData(fd) <- scramble, 'wellKey')

    expect_error(colData(fd) <- scramble[-1:-10,])
})

context('Testing data.table method')

test_that('Can cast to data.table', {
    dt <- as(fd, 'data.table')
    expect_is(dt, 'data.table')
    expect_equal(dt$value, as.vector(assay(fd)))
})

context('Play nicely with reshape/reshape2/data.table')
datArray <- array(c(1:39, NA), dim=c(2,4,5))
datList <- list(A=datArray, B=1:10)
test_that('Can melt with reshape2', {
    #try(detach('package:reshape', force=TRUE), silent=TRUE)
    tryCatch(library(reshape2, pos=length(search())), error = function(e) skip('Install reshape2'))
    M <- reshape2::melt(datArray, na.rm=TRUE, value.name='foo')
    expect_equal(M, melt(datArray, na.rm=TRUE, value.name='foo'))
    M2 <- reshape2::melt(datList, na.rm=TRUE, value.name='foo')
    expect_equal(M2, melt(datList, na.rm=TRUE, value.name='foo'))
    detach('package:reshape2')
})

test_that('Can melt with reshape', {
    library(reshape, pos=length(search()))
    M <- reshape::melt(datArray)
    expect_equal(M, melt(datArray))
    M2 <- reshape::melt(datList, value.name='foo')
    expect_equal(M2, melt(datList, value.name='foo'))
})

test_that('Can melt with data.table', {
    #library(data.table, pos=length(search()))
    dt <- data.table(A=rep(LETTERS, 2), B=rnorm(52))
    dtm <- melt(dt, id.var='A')
    expect_is(dtm, 'data.table')
})
