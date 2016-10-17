fd.small <- fd[,seq(from=1, to=ncol(fd), by=11)]
context('can construct FluidigmAssay')
test_that('ncells set', {
    expect_true('ncells' %in% names(colData(fd.small)))
})

context('Testing filtering')
test_that('Can filter',{
ff <- filter(fd.small, groups='ncells')
expect_that(ff, is_a('FluidigmAssay'))
ff2 <- filter(fd.small, groups='ncells', filt_control=list(nOutlier=1, sigmaContinuous=3)) #can't set sigmaContinuous to zero because we can't construct and empty object
expect_true(ncol(ff2)<ncol(ff))
ff2 <- filter(fd.small, groups='ncells', filt_control=list(nOutlier=1, sigmaContinuous=3, filter=FALSE), apply_filter=TRUE)
expect_that(ff2, is_a('list'))
ff3 <- filter(fd.small, apply_filter=FALSE)
expect_that(ff3, is_a('data.frame'))
})

context('Testing construction of reordered SingleCellAssay')
test_that('Can construct reorderd Assay',{
data(vbeta)
vbeta <- computeEtFromCt(vbeta)
# okay when data is ordered by wells, primers
vbeta.sca <- FromFlatDF(vbeta, idvars = c("Subject.ID", "Chip.Number", "Well"),
                             primerid = "Gene", measurement = "Et", geneid = "Gene",
                             cellvars = c("Number.of.Cells", "Population"), phenovars = c("Stim.Condition",
                                                                                          "Time"), id = "vbeta all")
expect_that(vbeta.sca,is_a("SingleCellAssay"))
# fails after data is ordered by primers,wells
reordering<-eval(as.call(c(order,as.list(vbeta[,c("Gene","Subject.ID", "Chip.Number", "Well")]))))
vbeta.reordered=vbeta[reordering,]
vbeta.reord.sca <- FromFlatDF(vbeta.reordered, idvars = c("Subject.ID", "Chip.Number", "Well"),
                                   primerid = "Gene", measurement = "Et", geneid = "Gene",
                                   cellvars = c("Number.of.Cells", "Population"), phenovars = c("Stim.Condition",
                                                                                                "Time"), id = "vbeta all reordered")
expect_that(vbeta.reord.sca,is_a("SingleCellAssay"))
})
