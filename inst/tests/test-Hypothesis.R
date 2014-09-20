context('Testing validity')
trms <- c('XXX', 'XXX1', 'YYY')
ch <- c('XXX', 'XXX1')
test_that('can construct', {
    manual <- new('Hypothesis', .Data=ch)
    expect_is(manual, 'Hypothesis')
    fun <- Hypothesis(ch)
    ## Get same result using static constructor or new
    expect_equal(manual, fun)
    fun <- generateHypothesis(manual, trms[1:2])
    expect_is(fun, 'Hypothesis')
    expect_equivalent(fun@transformed, diag(2))
    
    fun2 <- generateHypothesis(fun, trms)
    expect_that(fun, not(equals(fun2)))
    
    obj <- new('CoefficientHypothesis', 'XXX1')
    expect_is(obj, 'CoefficientHypothesis')
})
