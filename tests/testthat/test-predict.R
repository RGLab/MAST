context('Test prediction/imputation')
test_that('Predicted values match naive estimates', {
    data(vbetaFA)
    vbeta_scramble = vbetaFA[sample(nrow(vbetaFA)), sample(ncol(vbetaFA))]
    zlmVbeta <- zlm(~ Stim.Condition, subset(vbeta_scramble, select = ncells==1)[1:10,])
    nx = matrix(c(1, 0, 0, 1), nrow = 2,  ncol = 2)
    dimnames(nx)[[2]] =  colnames(coef(zlmVbeta, 'D'))
    pred = predict(zlmVbeta, modelmatrix = nx)
    expect_equal(pred[,muC], as.vector(coef(zlmVbeta, 'C')))
    expect_equal(pred[,etaD], as.vector(coef(zlmVbeta, 'D')))
})