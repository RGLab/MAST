context('Testing LRT')
n <- 100
w.x <- rep(1, n)
w.y <- rep(1, n)
x <- rnorm(n, 10)
y <- rnorm(n, 8)
SingleCellAssay:::lrtest(w.x, w.y, x, y)

w.x <- rep(1, n)
w.y <- rep(0, n)
x <- rnorm(n, 10)
y <- integer(0)
SingleCellAssay:::lrtest(w.x, w.y, x, y)

source('common-fixtures.R')
fd.spl <- split(fd, 'Number.of.Cells')
lrout <- SingleCellAssay:::lrt(fd.spl[[1]], 'Subject.ID', returnall=FALSE)
lrout2 <- suppressMessages(suppressWarnings(zlm.SingleCellAssay(et ~ Subject.ID, fd.spl[[1]], hypothesis.matrix='Subject.ID', type='LRT')))
context('testing for equality between glm lrtest and two-sample')
test_that('LRT and zlm are equivalent', {
expect_equivalent(lrout$lrstat, lrout2[,2,3])
})
