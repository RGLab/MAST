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
lrout <- lrt(fd.spl[[1]], 'Subject.ID', returnall=FALSE)
#lrout2 <- zlm.SingleCellAssay(et ~ Subject.ID, fd.spl[[1]])
#context('testing for equality between glm lrtest and two-sample')
#expect_equal(lrout$lrstat, lrout2[[1]]$lrstat)
