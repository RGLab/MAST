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

fd.spl <- split(fd, 'Number.of.Cells')
context('testing for equality between glm lrtest and two-sample')
test_that('LRT and zlm are equivalent', {
    lrout <- lrt(fd.spl[[1]], 'Subject.ID', returnall=FALSE)
    lrout2 <- suppressMessages(suppressWarnings(zlm.SingleCellAssay(~ Subject.ID, fd.spl[[1]], hypothesis='Subject.ID', type='LRT')))
    smallDOF <- freq(fd.spl[[1]])<=2/nrow(fd.spl[[1]])
    ## we are more conservative about returning low-DOF fits for dichotomous
    browser()
    expect_equivalent(lrout$lrstat[!smallDOF], lrout2[!smallDOF,3,1])
    expect_true(all(lrout$lrstat[smallDOF]>=lrout2[smallDOF,3,1]))
})
