context('Testing conversion from MASTClassic format')
test_that('Can convert', {
    old <- readRDS('SingleCellAssay-Object-MASTClassic.RData')
    newObj <- convertMASTClassicToSingleCellAssay(old)
    expect_is(newObj, 'SingleCellAssay')
    expect_equal(nrow(newObj), 10)
    expect_equal(ncol(newObj), 100)
})
