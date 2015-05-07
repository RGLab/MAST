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
    expect_equivalent(fun@contrastMatrix, diag(2))
    
    fun2 <- generateHypothesis(fun, trms)
    expect_that(fun, not(equals(fun2)))    
})

test_that('Coefficient Hypothesis Behaves',{
    obj <- new('CoefficientHypothesis', trms[1:2])
    expect_is(obj, 'CoefficientHypothesis')
    obj <- generateHypothesis(obj, trms[1:2])
    expect_equal(obj@index, 1:2)
    expect_equal(obj@contrastMatrix, diag(2), check.attributes=FALSE)
})

test_that('Can give terms at construction time',{
    delayed <- Hypothesis(ch)
    delayed <- generateHypothesis(delayed, trms)
    run <- Hypothesis(ch, trms)
    expect_equal(delayed, run)
})
