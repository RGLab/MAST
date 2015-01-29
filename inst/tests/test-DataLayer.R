dl0 <- new('DataLayer', array(0, dim=c(0, 75, 1)))

dl1 <- new('DataLayer', array(1, c(1, 1, 1), dimnames=list('A'=1,'B'=2, 'C'=3)))
dl <- new('DataLayer', array(1:16, c(2, 4, 2)))
dl2 <- new('DataLayer', array(1:8, c(2, 4, 1)))
mat <- matrix(1:8, nrow=2)

test_that('size zero dl preserves dimensions',{
  expect_equal(ncol(dl0), 75)
  expect_equal(nlayer(dl0), 1)
  expect_equal(nrow(dl0), 0)
})

context('nrow, ncol')

test_that('one by one dl works',{
  expect_equal(nrow(dl1), 1)
    expect_equal(ncol(dl1), 1)
  expect_false(is.null(dimnames(exprs(dl1))))
  expect_equal(dim(exprs(dl1)), c(1, 1))
})


test_that('singleton dl works', {
  expect_equal(nrow(dl2), 2)
  expect_equal(ncol(dl2), 4)
  expect_equal(nlayer(dl2), 1)
})

test_that('doubleton dl works', {
  expect_equal(exprs(dl), exprs(dl2))
  layer(dl) <- 2
  expect_false(all(exprs(dl)==exprs(dl2)))
  layer(dl) <- 1
})

context('subset and replace')
test_that('subset works', {
  expect_equal(dl[[,1]], mat[,1,drop=FALSE], check.names=FALSE)
  expect_equal(dl[[1,]], mat[1,,drop=FALSE], check.names=FALSE)
  ## repeated indices work as expected
  expect_equal(dl[[c(2, 2, 1, 2),]], mat[c(2, 2, 1, 2),,drop=FALSE], check.names=FALSE)
   expect_equal(dl[[,c(2, 2, 1, 2)]], mat[,c(2, 2, 1, 2),drop=FALSE], check.names=FALSE)
  matidx <- cbind(c(1, 1, 2), 1:3)
  expect_equal(dl[[matidx]], mat[matidx])
  mat[matidx] <- -999
  dl[[matidx]] <- -999
  expect_equal(exprs(dl), mat, check.attributes=FALSE)
})

test_that('exprs replace throws error for non-conforming', {
    expect_error(exprs(dl2) <- t(exprs(dl2)))
})

test_that('[ subscripting works', {
  expect_is(dl[1,], 'DataLayer')
  expect_is(dl[,2:3], 'DataLayer')
  expect_equal(exprs(dl[1,]), dl[[1,,drop=FALSE]], check.attributes=FALSE)
  expect_equal(exprs(dl[c(2, 2, 1, 2),]), dl[[c(2, 2, 1, 2),]], check.attributes=FALSE)
  expect_equal(exprs(dl[,c(2, 2, 1, 2)]), dl[[,c(2, 2, 1, 2)]], check.attributes=FALSE)
})

 spl <- split(dl, 1:2)
test_that('splitting works', {
 expect_is(spl, 'list')
 expect_equivalent(exprs(spl[[1]]), dl[[1,,drop=FALSE]])
})

test_that('combine works', {
  c <- combine(spl[[1]], spl[[2]])
  expect_is(c, 'DataLayer')
  expect_equivalent(c, dl)
})

context("Add layer")
test_that('Can add layer', {
    dl2 <- addlayer(dl, 'new')
    expect_true(all(is.na(exprsLayer(dl2, 'new'))))
    expect_equivalent(exprsLayer(dl2, 1), exprs(dl))

    exprsLayer(dl2, 'new') <- 999
    expect_true(all(exprsLayer(dl2, 'new')==999))
    layer(dl2) <- 'new'
    expect_true(all(exprs(dl2)==999))

})
