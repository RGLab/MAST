data(maits,package='MAST', envir = environment())
context('thresholding')
sca <- FromMatrix(t(maits$expressionmat[,1:1000]), maits$cdat, maits$fdat[1:1000,])

test_that('Can threshold', {
    tt <- thresholdSCRNACountMatrix(assay(sca))
    cp <-
        structure(c(2.55293247034713, 2.55293247034713, 4.89998792291031, 
                    4.83241868945595, 5.82780061591199, 2.25819983902936, 3.31358752487074
                    ), .Names = c("(0.0426,0.354]", "(0.354,0.757]", "(0.757,1.28]", 
                                  "(1.28,1.96]", "(1.96,2.84]", "(2.84,3.99]", "(3.99,13.2]"))
    expect_equal(cp, tt$cutpoint)
})

test_that('Warn when double-log transforming', {
    expect_warning(tt <- thresholdSCRNACountMatrix(assay(sca), data_log=FALSE), 'log')
})

test_that('Warn when erroneously not log-transforming', {
    expect_warning(try(thresholdSCRNACountMatrix(2^assay(sca)-1, data_log=TRUE) # will throw an error about no bimodal bins, too
                     , silent=TRUE), 'log')
    expect_error(hushWarning(thresholdSCRNACountMatrix(2^assay(sca)-1, data_log=TRUE), 'log'), 'bimodal')
})
