dl0 <- new('DataLayer', array(1, c(1, 1, 1), dimnames=list('A'=1,'B'=2, 'C'=3)))
dl <- new('DataLayer', array(1:16, c(2, 4, 2)))
dl2 <- new('DataLayer', array(1:8, c(2, 4, 1)))
mat <- matrix(1:8, nrow=2)
context('nrow, ncol')

test_that('one by one dl works',{
  expect_equal(nrow(dl0), 1)
    expect_equal(ncol(dl0), 1)
  expect_false(is.null(dimnames(exprs(dl0))))
  expect_equal(dim(exprs(dl0)), c(1, 1))
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
  expect_equal(dl[1,], mat[1,])
  matidx <- cbind(c(1, 1, 2), 1:3)
  expect_equal(dl[matidx], mat[matidx])
  mat[matidx] <- -999
  dl[matidx] <- -999
  expect_equal(exprs(dl), mat)
})
