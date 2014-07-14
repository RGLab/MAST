obj <- new('LMERlike', design=cData(fd), formula=~Stim.Condition + (1|Subject.ID))


if(require(lme4)){
    obj <- fit(obj, response=exprs(fd)[,2])
objC <- lmer(obj@response ~Stim.Condition +  (1|Subject.ID), data=obj@design, subset=obj@response>0, REML=FALSE)
objD <- glmer(obj@response>0 ~Stim.Condition + (1|Subject.ID), data=obj@design, family=binomial())
source('common-lmWrapper-tests.R', local=TRUE)
    detach('package:lme4')
}

