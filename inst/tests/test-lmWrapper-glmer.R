obj <- new('LMERlike', design=cData(fd), formula=~Stim.Condition + ncells + (1|Subject.ID))
if(require(lme4)){
source('common-lmWrapper-tests.R', local=TRUE)
}
