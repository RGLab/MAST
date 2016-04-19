library(testthat)

geneid="Gene"
primerid='Gene'
measurement='et'
idvars=c('Subject.ID', 'Chip.Number', 'Stim.Condition', 'Population', 'Well', 'Number.of.Cells')
phenovars=NULL
cellvars='Experiment.Number'
featurevars=NULL
ncells <- 'Number.of.Cells'

##Tests depending on vbeta
data(vbeta)
test_that("vbeta can be loaded",{
  expect_that(vbeta,is_a("data.frame"))
})

vbeta$et <- ifelse(is.na(vbeta$Ct), 0, 40-vbeta$Ct)


fd <- FromFlatDF(vbeta, idvars=idvars, primerid=primerid, measurement=measurement,cellvars=cellvars, geneid=geneid, ncells='Number.of.Cells', class='FluidigmAssay')
test_that('could create SingleCellAssay', {
    expect_that(fd, is_a('SingleCellAssay'))
})
