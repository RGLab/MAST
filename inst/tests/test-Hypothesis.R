context('Testing validity')
trms <- c('XXX', 'XXX1', 'YYY')
ch <- c('XXX=0', 'XXX1=0')
assign <- c(1, 1, 2)
test_that('can construct', {
    obj <- new('Hypothesis', characterHypothesis=ch, Terms=trms)
    expect_is(obj, 'Hypothesis')
    obj <- new('SimpleHypothesis', characterHypothesis=ch, Terms=trms)
    expect_is(obj, 'SimpleHypothesis')
    obj <- new('TermHypothesis', characterHypothesis=ch, Terms=trms, Assign=assign)
    expect_is(obj, 'TermHypothesis')
})


ch2 <- c('XXX=YYY')
ch3 <- c('XXX=0')
test_that('throw error on too general of hypothesis', {
    expect_error(new('SimpleHypothesis', characterHypothesis=ch2, Terms=trms), 'Terms must be individually or jointly compared to zero')
    expect_error(new('TermHypothesis', characterHypothesis=ch2, Terms=trms, Assign=assign), 'Terms must be individually or jointly compared to zero')
    expect_error(new('TermHypothesis', characterHypothesis=ch3, Terms=trms, Assign=assign), 'Must test all levels of a factor jointly')
})


trms2 <- c('XXX')
context('Accessors')
obj <- new('Hypothesis', characterHypothesis=ch, Terms=trms)
obj2 <- new('Hypothesis', characterHypothesis=ch3, Terms=trms2)
test_that('Can get contrast matrix', {
    expect_equal(dim(contrastMatrix(obj)), c(2, 3))
    expect_equal(dim(contrastMatrix(obj2)), c(1, 1))
})

test_that('Can get rhs', {
    expect_equal(dim(rhs(obj)), c(2, 1))
    expect_equal(dim(rhs(obj2)), c(1, 1))
})


context('Delayed construction')
test_that('Can delay construct', {
     obj <- new('Hypothesis', characterHypothesis=ch, Terms=trms, Assign=assign)
    expect_equivalent(obj,Hypothesis(ch)(trms, assign)[[1]])

     obj <- new('SimpleHypothesis', characterHypothesis=ch, Terms=trms, Assign=assign)
    expect_equivalent(obj,SimpleHypothesis(ch)(trms, assign)[[1]])

     obj <- new('TermHypothesis', characterHypothesis=ch, Terms=trms, Assign=assign)
    expect_equivalent(obj,TermHypothesis(ch)(trms, assign)[[1]])
     
})
