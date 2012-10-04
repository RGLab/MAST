doubleid <- data.frame(id1=1:3, id2=1:3, et=rep(3, 3), f1=rep('A', 3))
smallsc <- SingleCellAssay(doubleid, idvars='id1', geneid='f1', primerid='f1', measurement='et')

context('SingleCellAssay methods')
test_that('combine works', {
spl <- split(smallsc, 'id1')
c1 <- combine(spl[[1]], spl[[2]])
expect_that(c1, is_a('SingleCellAssay'))
expect_equal(nrow(c1), 2)
c2 <- combine(spl[[1]], spl[[2]], spl[[3]])
expect_that(c2, is_a('SingleCellAssay'))

c2 <- do.call(combine, spl@set)
expect_that(c2, is_a('SingleCellAssay'))
expect_equal(c2, smallsc)
})
