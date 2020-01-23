obj <- new('LMERlike', design=colData(fd), formula=~Stim.Condition + (1|Subject.ID), strictConvergence = FALSE)

context('LMERlike')
if(require(lme4)){
    obj <- fit(obj, response=t(assay(fd))[,2])
objC <- lmer(obj@response ~Stim.Condition +  (1|Subject.ID), data=as.data.frame(obj@design), subset=obj@response>0, REML=FALSE)
objD <- glmer(obj@response>0 ~Stim.Condition + (1|Subject.ID), data=as.data.frame(obj@design), family=binomial())
source('common-lmWrapper-tests.R', local=TRUE)

test_that('Signal error if no random effects', {
    data(vbetaFA)
    expect_error(zlm(~ Stim.Condition, vbetaFA[1:10,], method='glmer', ebayes = FALSE), 'specify at least one random effect') 
})

    test_that('lrt is non-NA', {
        z1 = zlm(~ Population + (1|Subject.ID), vbetaFA[3:5,], method='glmer', ebayes = FALSE)
        ## z2 = zlm(~ Population, vbetaFA[3:5,], method='bayesglm', ebayes = FALSE)
        
        ## VbetaResponsive is completely aliased with Subject.ID so test statistics should be no larger than the model ignoring subjectid
        ## And will be equal if the random effect is singular
        l1 = lrTest(z1, CoefficientHypothesis('PopulationVbetaResponsive'))
        expect_true(all(l1[,'hurdle', 'lambda']>0))
    })
    
    try(detach('package:lme4'), silent=TRUE)
}

