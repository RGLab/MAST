library(testthat)
library(devtools)
## Delete these lines before packaging
## library(flowCore)
## library(reshape)
## source('./R/AllClasses.R')
## source('R/SingleCellAssay-methods.R')
## source('R/SCASet-methods.R')
## source('R/Readers.R')
## End devel code

Files<-list.files(path=system.file("extdata/",package="SingleCellAssay",mustWork=TRUE),pattern="^sc.*csv$",full=T)[1]
test_that("Can find data files",{
  expect_that(Files,is_a("character"))
})
dat<- read.csv(Files, as.is=TRUE)
idvars = c("Patient.ID", "Stim.Agent", "Chip.Number", "Well")
dat$ct <- ifelse(is.na(dat$X40.Ct), 0, dat$X40.Ct)


dat_complete <- subset(dat, Patient.ID != '200-143')

geneid="Assay.Name"
primerid='Assay.Name'
measurement='ct'

sc <- SingleCellAssay(dat_complete, idvars=idvars, geneid=geneid, primerid=primerid, measurement=measurement)

themap <- new("Mapping",mapping=list(idvars=idvars, geneid=geneid, primerid=primerid, measurement=measurement))
cellvars <- c('Date.of.Sort', 'Sero.Status', 'Time.of.Stim', 'Stim.Agent', 'Number.of.Cells')
featurevars <- c('Assay.Catalogue', 'Assay.Module')
phenovars <- NULL


##Tests depending on vbeta
data(vbeta)
test_that("vbeta can be loaded",{
  expect_that(vbeta,is_a("data.frame"))
})

