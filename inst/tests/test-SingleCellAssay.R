source('common-fixtures.R') #gets us Files, dat, dat_complete, sc

context("Loading Mario's Data")
test_that("csv could be read",{
  expect_that(dat,is_a("data.frame"))
  expect_that( all(c(idvars, 'ct') %in% names(dat)), is_true())
})


context("Generating a complete and incomplete subset")
countComplete <- table(do.call(paste, dat_complete[,idvars]))
expect_that(all(countComplete==countComplete[1]), is_true())
dat_incomplete <- dat_complete[-seq(1, nrow(dat_complete), by=101),]
counts <- table(do.call(paste, dat_incomplete[,idvars]))
expect_that(all(counts == counts[1]), is_false())

blank <- dat_complete[1,]
blank$`__wellKey` <- 1
blankinst <- SingleCellAssay(blank, idvars=idvars, geneid="Assay.Name", primerid='Assay.Name', measurement='ct')
test_that("Can create empty instance with pre-existing wellKey",{
  expect_that(blankinst, is_a("SingleCellAssay"))
})

test_that("Can get wellKey",{
  expect_that(getwellKey(blankinst)[[1]], equals(digest(paste(blankinst@env$data[,getMapping(blankinst@mapping)$idvars],collapse=" "))))
})
          
test_that("Can load complete data", {
  expect_that(sc, is_a("SingleCellAssay"))
  tab <- table(melt(sc)$`__wellKey`)
  expect_that(tab, is_equivalent_to(countComplete))
})

test_that("Cellkey unique identifies a cell", {
  tab <- table(melt(sc)$`__wellKey`, do.call(paste, melt(sc)[, idvars]))
  expect_true(all(tab %in% c(0,96)))
  
})

sci<- SingleCellAssay(dat_incomplete, idvars=idvars, geneid="Assay.Name", primerid='Assay.Name', measurement='ct')
test_that("Completes incomplete data", {
  expect_that(sci, is_a("SingleCellAssay"))
  expect_equal(nrow(melt(sci)), nrow(dat_complete))
})

context('testing cell and feature dictionaries')

scd <- SingleCellAssay(dat_complete, mapping=themap, featurevars=featurevars, cellvars=cellvars, phenovars=phenovars)

test_that('Cell data has correct number of row/columns', {
  expect_that(cData(scd), is_a('data.frame'))
  expect_equivalent(nrow(scd), nrow(cData(scd)))
  expect_equal(unique(cData(scd)), cData(scd))
  expect_equivalent(ncol(cData(scd)), length(unique(c(phenovars, idvars, cellvars))))
})

test_that('Feature data has correct number of row/columns', {
 expect_that(fData(scd), is_a('data.frame'))
  expect_that(featureData(scd), is_a('AnnotatedDataFrame'))
   expect_equivalent(ncol(scd), nrow(fData(scd)))
  expect_equal(unique(fData(scd)), fData(scd))
  expect_equivalent(ncol(fData(scd)), length(unique(c(featurevars, geneid, primerid))) )
})

context("Testing methods")
test_that("Can subset complete data with integer indices",{
  ind<- seq(1, length(countComplete), by=5) 
  sub <- sc[[ind]]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_that(getwellKey(sub), equals(getwellKey(sc)[ind]))
})

  boolind<- rep(c(FALSE, TRUE, FALSE), length.out=length(countComplete))
test_that("Can subset complete data with boolean indices",{
  sub <- sc[[boolind]]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_that(getwellKey(sub), equals(getwellKey(sc)[which(boolind)]))
})

test_that('Subsetting preserves cell and featuredata',{
  sub <- scd[[boolind]]
  expect_that(sub, is_a("SingleCellAssay"))
  expect_that(getwellKey(sub), equals(getwellKey(scd)[which(boolind)]))
  expect_equal(featureData(sub), featureData(scd))
  expect_equivalent(cData(sub), cData(scd)[boolind,])
  sub2 <- sub[[boolind]]
  expect_equivalent(cData(sub2), cData(sub)[boolind,])
})

exprComplete <- exprs(sc)
test_that("Exprs works", {
  expect_is(exprComplete, "matrix")
  expect_equal(nrow(exprComplete), length(getwellKey(sc)))
  ind <- seq(1, nrow(dat_complete), by=1042)
  expect_equal(melt(sc)$ct[ind], exprComplete[ind])
  geneandrow <- melt(sc)[1054,c("Assay.Name", "__wellKey")]  
  thect <- melt(sc)[1054, "ct"]
  expect_equal(exprComplete[geneandrow[[2]], geneandrow[[1]]], thect)
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
  stim <- table(cData(scd)$Stim.Agent)[1]
  sub1 <- subset(scd, Stim.Agent == names(stim))
  expect_equivalent(nrow(sub1), stim)
})

test_that('Subset throws an intelligent error if thesubset cannot be evaluated', {
 expect_that(subset(scd, NOTPRESENT==fdsfjkl), throws_error('not found'))
})

context("SCASet works")
test_that('Can construct', {
aset <- SCASet(melt(scd), splitby='Patient.ID', mapping=getMapping(scd))
})

test_that('Can split',{
  splat <- split(scd, cData(scd)$Patient.ID)
  expect_that(splat, is_a('SCASet'))
  splat.byfieldname <- split(scd, 'Patient.ID')
    expect_that(splat.byfieldname, is_a('SCASet'))
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

# test_that("Can add mappings",{
#   addMapping(sc, c("plate"="Chip.Number"))
# })
# 
# test_that("Can rewrite mappings", {
#   addMapping(sc, c("plate" = "Date.of.Sort"))
# })


## context('testing error handling')
## badmap <- themap
## badmap$idvars <- c(badmap$idvars, 'Assay.Name')
## test_that('We throw an error when idvars is not disjoint from probeid',{
##   expect_that(SingleCellAssay(dat_complete, mapping=badmap), throws_error())
## })
