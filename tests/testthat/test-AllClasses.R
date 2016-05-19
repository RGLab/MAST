singleid <- data.frame(id1=1:3, et=rep(3, 3), f1=rep('A', 3))
doubleid <- data.frame(id1=1:3, id2=LETTERS[1:3], et=rep(3, 3), f1=rep('A', 3))
duplicateprimers <- data.frame(id1=rep(1, 5), et=rep(3, 5), f1=c(rep('B', 1), rep('A', 2), rep('C', 2)))
duplicateprimers2 <- data.frame(id1=rep(1, 5), et=rep(3, 5), primerid=c(rep('B', 1), rep('A', 2), rep('C', 2)))
integerprimers <- data.frame(id1=rep(1, 3), et=rep(3, 3), f1=1:3)


context('Testing if we use drop=FALSE correctly')
test_that('Create a two-column id SingleCellAssay', {
  tmp <- FromFlatDF(doubleid, idvars=c('id1', 'id2'), geneid='f1', primerid='f1', measurement='et')
  expect_that(tmp, is_a('SingleCellAssay'))
})

intprimer <- FromFlatDF(integerprimers, idvars='id1', geneid='f1', primerid='f1', measurement='et')
test_that('Create  SingleCellAssay from integerprimers', {
  expect_that(intprimer, is_a('SingleCellAssay'))
  expect_is(mcols(intprimer)$primerid, 'character')
})


tmp <- FromFlatDF(singleid, idvars='id1', geneid='f1', primerid='f1', measurement='et')
test_that('Create a one-column id SingleCellAssay', {
    expect_that(tmp, is_a('SingleCellAssay'))
})

test_that('duplicate primers doesn\'t throw error', {
  expect_warning(sc <- FromFlatDF(dataframe=duplicateprimers, idvars='id1', geneid='f1', primerid='f1', measurement='et'), 'A')
  expect_equal(nrow(sc), 5)
  expect_warning(sc <- FromFlatDF(dataframe=duplicateprimers2, idvars='id1', primerid='primerid', measurement='et'), 'A')
  expect_true('primerid.orig' %in% names(mcols(sc)))
})

