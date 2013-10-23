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

context('Can construct')
cdat <- as(data.frame(wellKey=c('A', 'B'), ncells=c(1, 1)), 'AnnotatedDataFrame')
fdat <- as(data.frame(primerid=LETTERS[1:4], meta=TRUE), 'AnnotatedDataFrame')
val <- new('DataLayer', array(1:8, c(2, 4, 1), dimnames=list(wellKey=cdat$wellKey, primerid=fdat$primerid, measurement='et')))
sc <- new('SingleCellAssay', .Data=val, cellData=cdat, featureData=fdat)
fd <- new('FluidigmAssay', .Data=val, cellData=cdat, featureData=fdat)


context("Generating a complete and incomplete subset")
dat_complete <- vbeta
countComplete <- table(do.call(paste, dat_complete[,idvars]))
expect_that(all(countComplete==countComplete[1]), is_true())
dat_incomplete <- dat_complete[-seq(1, nrow(dat_complete), by=101),]
counts <- table(do.call(paste, dat_incomplete[,idvars]))
expect_that(all(counts == counts[1]), is_false())

blank <- dat_complete[1,]
blankinst <- new('SingleCellAssay', dataframe=blank, idvars=idvars, primerid=geneid, measurement=measurement)
blankinst.fd <- FluidigmAssay(dataframe=blank, idvars=idvars, primerid=geneid, measurement=measurement, ncells=ncells)
test_that("Can create empty instance with pre-existing wellKey",{
  expect_that(blankinst, is_a("SingleCellAssay"))
})

## test_that("Can get wellKey",{
##   expect_that(getwellKey(blankinst)[[1]], equals(digest(paste(blankinst@env$data[,getMapping(blankinst@mapping)$idvars],collapse=" "))))
## })

vbeta$et <- ifelse(is.na(vbeta$Ct), 0, 40-vbeta$Ct)
fd <- FluidigmAssay(vbeta, idvars=idvars, primerid=primerid, measurement=measurement, ncells=ncells, geneid=geneid, keep.names=TRUE)
test_that('could create FluidigmAssay', {
  expect_that(fd, is_a('SingleCellAssay'))
    expect_that(fd, is_a('FluidigmAssay'))
})


context('Testing legacy behavior')
test_that('sort unsorted SingleCellAssay', {
  ss <- fd[,c(5, 1, 4, 10, 15)]         #or should we return a sorted object here?
  sorted <- new('SingleCellAssay', .Data=as(ss, 'DataLayer'), cellData=cellData(ss), featureData=featureData(ss), sort=TRUE)
  expect_equal(sort(fData(sorted)$primerid), fData(sorted)$primerid)
})
test_that('keep field names',{
  m <- melt(fd)
  expect_true(all(c(primerid, measurement, ncells, geneid) %in% names(m)))
})



sc <- fd
test_that("Can load complete data", {
  expect_that(sc, is_a("SingleCellAssay"))
  tab <- table(melt(sc)$wellKey)
  expect_that(tab, is_equivalent_to(countComplete))
})

test_that("Cellkey unique identifies a cell", {
  tab <- table(SingleCellAssay:::melt(sc)$wellKey, do.call(paste, SingleCellAssay:::melt(sc)[, idvars]))
  expect_true(all(tab %in% c(0,75)))
  
})


context('test construction helper funcs')
test_that("uniqueModNA doesn't include NA", {
  naframe <- data.frame(var=rep(c(1, 2), each=2), na=c(NA, -9, NA, -9))
  expect_equal(nrow(uniqueModNA(naframe, exclude='var')), 2)
  expect_equal(nrow(as.data.frame(uniqueModNA(naframe[,-2, drop=FALSE], exclude='var'))), 2)
})

sci<- SingleCellAssay(dat_incomplete, idvars=idvars, primerid=geneid, measurement=measurement)
test_that("Completes incomplete data", {
  expect_that(sci, is_a("SingleCellAssay"))
  expect_equal(nrow(SingleCellAssay:::melt(sci)), nrow(dat_complete))

  incomplete <- rbind(melt(fd[1:20,1:20]),
                      melt(fd[21:50, 11:30])) #equally sized primerid blocks
  fd.incomplete <- FluidigmAssay(incomplete, idvars=idvars, primerid=primerid, measurement=measurement, ncells='ncells', geneid=geneid, keep.names=TRUE)
  expect_message(FluidigmAssay(incomplete, idvars=idvars, primerid=primerid, measurement=measurement, ncells='ncells', geneid=geneid, keep.names=TRUE), 'incomplete')
  expect_equal(nrow(fd.incomplete), 50)
  expect_equal(ncol(fd.incomplete), 30)
  expect_true(any(is.na(exprs(fd.incomplete))))
  expect_equal(exprs(fd.incomplete[21:50, 11:30]),exprs(fd[21:50, 11:30]))
  expect_equal(exprs(fd.incomplete[1:20, 1:20]),exprs(fd[1:20, 1:20]))
})

## No more mapping, hurray!
## context('testing cell and feature dictionaries')

## complete<-(SingleCellAssay:::melt(sc))[,setdiff(colnames(SingleCellAssay:::melt(sc)),"__wellKey")]
## scd <- SingleCellAssay(complete, mapping=getMapping(sc))

## test_that('Cell data has correct number of row/columns', {
##   expect_that(cData(scd), is_a('data.frame'))
##   expect_equivalent(nrow(scd), nrow(cData(scd)))
##   expect_equal(unique(cData(scd)), cData(scd))
##   expect_equivalent(ncol(cData(scd)), length(unique(c(phenovars, idvars, cellvars))))
## })
scd <- new('SingleCellAssay', .Data=sc@.Data, cellData=cellData(sc), featureData=featureData(sc))
test_that('Feature data has correct number of row/columns', {
 expect_that(fData(scd), is_a('data.frame'))
  expect_that(featureData(scd), is_a('AnnotatedDataFrame'))
   expect_equivalent(ncol(scd), nrow(fData(scd)))
  expect_equal(unique(fData(scd)), fData(scd))
 #one extra column for the primerid
  expect_equivalent(ncol(fData(scd)), length(unique(c(featurevars, geneid, primerid)))+1) 
})

context("Testing methods")
test_that("Can subset complete data with integer indices",{
  ind<- seq(1, length(countComplete), by=5) 
  sub <- sc[[ind,]]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_that(getwellKey(sub), equals(getwellKey(sc)[ind]))
  sub <- sc[[ind,"B3GAT1"]]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_that(getwellKey(sub), equals(getwellKey(sc)[ind]))
})

  boolind<- c(FALSE, TRUE, FALSE)
test_that("Can subset complete data with boolean indices",{
  sub <- sc[[boolind]]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_that(getwellKey(sub), equals(getwellKey(sc)[boolind]))
})

## This makes more sense to me than to propagate new wells/features consisting entirely of NA
test_that('NAs throw error when subsetting', {
  expect_error(sc[, c(1:4, NA)])
  expect_error(sc[c(boolind, NA), ])
})

test_that("Throw error when indexing with factors", {
    expect_error(sc[, factor(c('B3GAT1'))])
    expect_error(sc[factor('A'),])
})

test_that("loading Matrix package doesn't clobber generic table", {
  ## Matrix dispatches [[ and [ on x, i, j, drop
  ## The distinction between ANY and missing doesn't matter for
  ## i and j until Matrix is loaded due to method caching
  ## but then the dispatch occurs separately and things can break
  library(Matrix)
  sub <- sc[[boolind]]
  sub2 <- sc[boolind]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_equal(sub, sub2)
  detach('package:Matrix')
})

test_that('can subset by character', {
  sub <- sc[,'GAPDH']
  expect_equal(ncol(sub), 1)
  expect_error(sc[,'NOT_PRESENT'], 'NOT_PRESENT')
  expect_error(sc['NOT_PRESENT',], 'NOT_PRESENT')
  expect_equal(getwellKey(fd["Sub02 3 Stim(SEB) VbetaResponsive C02",]), "Sub02 3 Stim(SEB) VbetaResponsive C02")
})

test_that('Subsetting preserves cell and featuredata',{
  sub <- sc[[boolind]]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_that(getwellKey(sub), equals(getwellKey(scd)[boolind]))
  expect_equal(featureData(sub), featureData(scd))
  expect_equivalent(cData(sub), cData(scd)[boolind,])
#expect_equivalent(cData(sub2), cData(sub)[boolind,]) #check if the cData subset is the same as the subset of the cData
})

exprComplete <- exprs(sc)
test_that("Exprs works", {
  measurement <- 'et'                #fix so melt renames column correctly
  expect_is(exprComplete, "matrix")
  expect_equal(nrow(exprComplete), length(getwellKey(sc)))
  ind <- seq(1, nrow(dat_complete), by=1042)
  expect_equal(melt(sc)[ind,measurement], as.vector(exprComplete)[ind])
  geneandrow <- melt(sc)[1054,c(geneid, "wellKey")]  
  thect <- melt(sc)[1054, measurement]
  expect_equivalent(exprComplete[geneandrow[[2]], geneandrow[[1]]], thect)
})

test_that('Subset with TRUE is unchanged', {
  suball <- subset(scd, TRUE)
  expect_equal(suball, scd)
})

test_that('Subset with FALSE returns empty set', {
  subnone <- subset(scd, FALSE)
  expect_that(all.equal(scd, subnone), is_a('character'))
  expect_equal(nrow(subnone), 0)
})


test_that('Subset with names from SingleCellAssay works', {
  stim <- table(cData(scd)$Stim.Condition)[1]
  sub1 <- subset(scd, Stim.Condition == names(stim))
  expect_equivalent(nrow(sub1), stim)
})

test_that('Subset throws an intelligent error if thesubset cannot be evaluated', {
 expect_that(subset(scd, NOTPRESENT==fdsfjkl), throws_error('not found'))
})

context("SCASet works")
test_that('Can construct', {
aset <- SCASet(melt(scd), primerid=primerid, idvars=idvars, measurement='et', splitby='Subject.ID')
})

splat <- split(scd, cData(scd)$Subject.ID)
test_that('Can split',{
  expect_that(splat, is_a('SCASet'))
  splat.byfieldname <- split(scd, 'Subject.ID')
    expect_that(splat.byfieldname, is_a('SCASet'))
    splat <- split(scd, c('Subject.ID', 'Population'))
  expect_that(splat, is_a('SCASet'))
    splat <- split(scd, list(factor(cData(scd)$Subject.ID), factor(cData(scd)$Population)))
  expect_that(splat, is_a('SCASet'))
  
})

test_that('Can coerce to/from list', {
   tolist <- as(splat, 'list')
  expect_that(tolist, is_a('list'))
  expect_that(as(tolist, 'SCASet'), is_a('SCASet'))
})

context('Copy and replace')
test_that('Replace works', {
 exprs(scd) <- -111
 expect_true(all(exprs(scd)==(-111)))
})

test_that('Copy works', {
 sc2 <- copy(scd)
 expect_that(scd, equals(sc2))
 exprs(sc2) <- -999
 expect_false(any(exprs(scd) == -999))
})


context('Combine works')
doubleid <- data.frame(id1=1:3, id2=1:3, et=rep(3, 3), f1=rep('A', 3))
smallsc <- SingleCellAssay(doubleid, idvars='id1', geneid='f1', primerid='f1', measurement='et', id='1')

test_that('combine works', {
spl <- split(smallsc, 'id1')
c1 <- combine(spl[[1]], spl[[2]])
expect_that(c1, is_a('SingleCellAssay'))
expect_equal(nrow(c1), 2)
c2 <- combine(spl[[1]], spl[[2]], spl[[3]])
expect_that(c2, is_a('SingleCellAssay'))
expect_that(combine(spl), is_a('SingleCellAssay'))
})

test_that('combine throws error for non-conforming',{
  expect_error(combine(fd, spl[[1]]))
})






context('Testing FluidigmAssay')
fd <- as(scd, 'FluidigmAssay')
back <- as(fd, 'SingleCellAssay')
test_that('Can cast', {
  expect_is(fd, "FluidigmAssay")
  expect_is(back, 'SingleCellAssay')
  expect_equivalent(back, scd)
})

test_that('Can split FluidigmAssays', {
 splat.byfieldname <- split(fd, 'Subject.ID')
    expect_that(splat.byfieldname, is_a('SCASet'))
})

context('Test replace methods')


test_that('Can replace cData', {
    cDat <- cData(fd)
    cDat$foo <- "bar"
    cData(fd) <- cDat
    expect_true('foo' %in% names(cData(fd)))

    empty <- data.frame()
    expect_error(cData(fd) <- empty, 'wellKey')

    scramble <- cDat[sample(nrow(cDat)),]
    expect_warning(cData(fd) <- scramble, 'sorting')

    expect_error(cData(fd) <- scramble[-1:-10,], 'missing some wellkeys')
})
