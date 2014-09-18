test_that('Ebayes shrinkage corner cases', {
    ## zero DOF
    obj@priorDOF <- 0
    obj <- fit(obj)
    expect_equal(obj@fitC$dispersion, summary(objC)$dispersion)
    ## Complete shrinkage
    obj@priorDOF <- 1e20
    obj@priorVar <- 999
    obj <- fit(obj)
    expect_equal(obj@fitC$dispersion, 999, tol=1e-5)
    expect_equal(obj@fitC$dispersionNoShrink, summary(objC)$dispersion)
})

context('Torture test for convergence')
test_that('Throw error or signal non-convergence', {
    obj3 <- new(class(obj), design=data.frame(x=rep(1, 4)), formula=~1, fitArgsD=list(start=-1.81))
    ## An error is acceptable
     tt <- try({
    obj3 <- fit(obj3, response=c(1,1,1,0))
})
    if(!is(tt, 'try-error')){
        expect_true(abs(coef(obj3, 'D'))<2)
    }
})
