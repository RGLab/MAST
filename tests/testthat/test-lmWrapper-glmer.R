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
    try(detach('package:lme4'), silent=TRUE)
}

