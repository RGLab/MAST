data(maits,package='MAST', envir = environment())
context('thresholding')
sca <- FromMatrix(t(maits$expressionmat[,1:1000]), maits$cdat, maits$fdat[1:1000,])

test_that('Cutpoints reordered correctly', {
    somecp = c(4, 1:4, 3)
    expect_equal(orderCutpoints(1, somecp), rep(4, length(somecp)))
    expect_equal(orderCutpoints(length(somecp), somecp), c(1, 1, 2, 3, 3, 3))
    expect_equal(orderCutpoints(2, somecp), c(1, 1, 2, 3, 4, 4))
})

test_that('Can threshold', {
    tt <- thresholdSCRNACountMatrix(assay(sca))
    expect_equal(tt$cutpoint, sort(tt$cutpoint))
    cp <-
        structure(c(2.25819983902936, 2.25819983902936, 2.25819983902936, 2.25819983902936, 
                    2.25819983902936, 2.25819983902936, 3.31358752487074),
                    .Names = c("(0.0426,0.354]", "(0.354,0.757]", "(0.757,1.28]", 
                                  "(1.28,1.96]", "(1.96,2.84]", "(2.84,3.99]", "(3.99,13.2]"))
    expect_equal(cp, tt$cutpoint)
})

test_that('Warn when double-log transforming', {
    expect_warning(tt <- thresholdSCRNACountMatrix(assay(sca), data_log=FALSE), 'log')
})

test_that('Warn when erroneously not log-transforming', {
    expect_warning(try(thresholdSCRNACountMatrix(2^assay(sca)-1, data_log=TRUE) # will throw an error about no bimodal bins, too
                     , silent=TRUE), 'log')
})
