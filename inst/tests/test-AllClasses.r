singleid <- data.frame(id1=1:3, et=rep(3, 3), f1=rep('A', 3))
doubleid <- data.frame(id1=1:3, id2=1:3, et=rep(3, 3), f1=rep('A', 3))
duplicateprimers <- data.frame(id1=rep(1, 5), et=rep(3, 5), f1=c(rep('B', 1), rep('A', 2), rep('C', 2)))
duplicateprimers2 <- data.frame(id1=rep(1, 5), et=rep(3, 5), primerid=c(rep('B', 1), rep('A', 2), rep('C', 2)))

context('Testing if we use drop=FALSE correctly')
test_that('Create a two-column id SingleCellAssay', {
  tmp <- SingleCellAssay(doubleid, idvars=c('id1', 'id2'), geneid='f1', primerid='f1', measurement='et')
  expect_that(tmp, is_a('SingleCellAssay'))
})

tmp <- SingleCellAssay(singleid, idvars='id1', geneid='f1', primerid='f1', measurement='et')
test_that('Create a one-column id SingleCellAssay', {
    expect_that(tmp, is_a('SingleCellAssay'))
})

test_that('duplicate primers doesn\'t throw error', {
  expect_warning(sc <- SingleCellAssay(dataframe=duplicateprimers, idvars='id1', geneid='f1', primerid='f1', measurement='et'), 'A')
  expect_equal(ncol(sc), 5)
  expect_warning(sc <- SingleCellAssay(dataframe=duplicateprimers2, idvars='id1', primerid='primerid', measurement='et'), 'A')
  expect_true('primerid.orig' %in% names(fData(sc)))
})

context('Construction from matrix')
test_that('Can recreate', {
    tmp2 <- suppressMessages( FromMatrix('SingleCellAssay', exprs(tmp), cData(tmp), fData(tmp)))
    expect_is(tmp2, 'SingleCellAssay')
    expect_equivalent(tmp, tmp2)
    cd <- cData(tmp)
    cd$ncells <- 1
    expect_message(tmp3 <- FromMatrix('FluidigmAssay', exprs(tmp), cd), 'dimnames')
    expect_is(tmp3, 'FluidigmAssay')
})
