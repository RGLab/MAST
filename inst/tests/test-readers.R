source('common-fixtures.R')

context("Loading Raw Data")
test_that("testthat is loaded",{
expect_that(require(testthat),is_true())
#expect_that(require(SingleCellAssay),is_true())
expect_that(require(gdata),is_true())
expect_that(require(reshape),is_true())
})
F<-list.files(path=system.file("extdata/",package="SingleCellAssay",mustWork=TRUE),pattern="cmv.*.xls$",full=T)[1:5]
test_that("Can find data files",{
expect_that(F,is_a("character"))
})
	metadatafile<-list.files(path=system.file("extdata/",package="SingleCellAssay",mustWork=TRUE),pattern="metadata.csv$",full=T)
test_that("Can find metadata",{
expect_that(metadatafile,is_a("character"))
})
data<-read.fluidigm(F,metadata=metadatafile,metadataColClasses=c("character","factor","factor","factor","factor","factor"),meta.key="Filename",idvars=c("sample","panel","plate"),splitby="donor")
test_that("data can be loaded",{
expect_that(data,is_a("SCASet"))
})
test_that("SCASet can be subset with [[",{
expect_that(data[[1]],is_a("SingleCellAssay"))
})
test_that("SCASet works with lapply",{
expect_that(lapply(data,function(x) x)[[1]],is_a("SingleCellAssay"))
})
test_that("exprs works on SingleCellAssay",{
expect_that(exprs(data[[1]]),is_a("matrix"))
})


context("Loading Mario's Data")
F<-list.files(path=system.file("extdata/",package="SingleCellAssay",mustWork=TRUE),pattern="^sc.*csv$",full=T)[1:4]
test_that("Can find data files",{
expect_that(F,is_a("character"))
})
	metadatafile<-list.files(path=system.file("extdata/",package="SingleCellAssay",mustWork=TRUE),pattern="csv",full=T)
test_that("Can find metadata",{
expect_that(metadatafile,is_a("character"))
})
data<-read.fluidigm(files=F,metadata=NULL,meta.key=NULL,idvars=c("Sample.ID","Stim.Condition","Stim.Agent"),splitby="Patient.ID",unique.well.id="Well",raw=FALSE,gene="Assay.Name",measurement.processed="40-Ct")
test_that("data can be loaded",{
expect_that(data,is_a("SCASet"))
})
test_that("SCASet can be subset with [[",{
expect_that(data[[1]],is_a("SingleCellAssay"))
})
test_that("SCASet works with lapply",{
expect_that(lapply(data,function(x)x)[[1]],is_a("SingleCellAssay"))
})
test_that("exprs works on SingleCellAssay",{
expect_that(exprs(data[[1]]),is_a("matrix"))
})
