context('Construction from matrix')
#library(S4Vectors) #This shouldn't be necessary
NC <- 10
NF <- 100
assay <- matrix(rnorm(NC*NF), ncol=NC, nrow=NF)
assayDn <- matrix(rnorm(NC*NF), ncol=NC, nrow=NF,  dimnames=list(paste('p', 1:NF), LETTERS[1:NC]))
colData <- DataFrame(dat=runif(NC))
colDataRn <- DataFrame(colData, row.names=LETTERS[1:NC])
rData <- DataFrame(pid= 1:NF, row.names=paste('pid', 1:NF))
test_that('Can create', {
    se <- suppressMessages(FromMatrix(assay, colData, rData))
    expect_is(se, 'SummarizedExperiment0')
    expect_false(is.null(dimnames(se)[[1]]))
    expect_false(is.null(dimnames(se)[[2]]))
})

test_that('Preserve dimnames in exprsArray', {
    expect_silent(se2 <- FromMatrix(assayDn))
    expect_equal(dimnames(se2), dimnames(assayDn))
})

test_that('Throw error for mismatching names', {
    rData[['primerid']] <- 1:NF
    expect_error(FromMatrix(assayDn, colData, rData))
})

test_that('Integer primerids cast to character', {
    rData[['primerid']] <- 1:NF
    se <- FromMatrix(assay, colData, rData)
    expect_is(mcols(se)$primerid, 'character')
})
