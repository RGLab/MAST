context('Testing LRT')
n <- 100
w.x <- rep(1, n)
w.y <- rep(1, n)
x <- rnorm(n, 10)
y <- rnorm(n, 8)
MAST:::lrtest(w.x, w.y, x, y)

w.x <- rep(1, n)
w.y <- rep(0, n)
x <- rnorm(n, 10)
y <- integer(0)
MAST:::lrtest(w.x, w.y, x, y)

context('testing for equality between glm lrtest and two-sample')
test_that('LRT and zlm are equivalent', {
    fd.spl <- split(fd, 'Number.of.Cells')
    lrout <- lrt(fd.spl[[1]], 'Subject.ID', returnall=FALSE)
    hushWarning(zlm2 <- zlm.SingleCellAssay(~ Subject.ID, fd.spl[[1]], silent=FALSE), '(At least one component)|(No positive observations)')
    lrout2 <- lrTest(zlm2, 'Subject.ID')
    smallDOF <- freq(fd.spl[[1]])<=3/ncol(fd.spl[[1]]) | (1-freq(fd.spl[[1]]))<(3/ncol(fd.spl[[1]]))
    ##if(!isTRUE(all.equal(lrout$lrstat[!smallDOF], lrout2[!smallDOF,3,1], tolerance=1e-6, check.attributes=FALSE))) browser() 
    ## we are more conservative about returning low-DOF fits for dichotomous
    expect_equal(lrout$lrstat[!smallDOF], lrout2[!smallDOF,3,1], tol=1e-6, check.attributes=FALSE)
    expect_true(all(lrout$lrstat[smallDOF]>=lrout2[smallDOF,3,1]))
})
