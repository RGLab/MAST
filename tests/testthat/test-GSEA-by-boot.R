context('GSEA')
data(vbetaFA)
library(plyr)
vb1 = subset(vbetaFA[1:24,], ncells==1)
#vb1 = vb1[,freq(vb1)>.1]
zf = zlm(~Stim.Condition, vb1)
set.seed(1234)
boots = bootVcov1(zf, R = 36)
## replace NAs for each coefficient, gene and component
bootsUncor <- apply(boots, 2:4, function(col){
    col[!is.na(col)] <- col[!is.na(col)]-mean(col[!is.na(col)], na.rm=TRUE)
    col[is.na(col)] <- 0
    col
})

bootsUncor <- aaply(bootsUncor, 3:4, function(mat){
    S <- svd(mat, nv=ncol(mat))
    mat %*% S$v
})
bootsUncor <- aperm(bootsUncor, c(3,4,1,2))
sets=list(A=1:5, B=3:10, C=15, D=1:5, E=12:14, F=15:24, G=12, H=13, I=14)
sets[['notE']] <- setdiff(1:nrow(vb1), sets[['E']])
gseaClass <- gseaAfterBoot(zf, bootsUncor, sets, CoefficientHypothesis('Stim.ConditionUnstim'), control=list(n_randomize=Inf))
gsea <- gseaClass@tests
test_that('equal sets yield equal results', {
    expect_equal(gsea['A',,,],gsea['D',,,])
})

test_that('Singletons agree with coefficients', {
    expect_equal(gsea['C','cont','stat','test'], coef(zf, 'C')[15,'Stim.ConditionUnstim'])
    expect_equal(gsea['C','disc','stat','test'], coef(zf, 'D')[15,'Stim.ConditionUnstim'])
})

## 
test_that('Triples work as expected', {
              tripleEachT <- gsea[c('G', 'H', 'I'), 'disc', 'stat', 'test']
              tripleEachVar <- gsea[c('G', 'H', 'I'), 'disc', 'var', 'test']
              tripleAvgT <- sum(tripleEachT)/3
              expect_equal(tripleAvgT, gsea[c('E'), 'disc', 'stat', 'test'])
              tripleAvgVar <- sum(tripleEachVar)/9
              expect_equal(tripleAvgT, gsea[c('E'), 'disc', 'stat', 'test'])
          })

test_that('Avg cor is nearly null', {
    ## We are rank deficient or something when whitening the correlation for this set
    ## in the continuous component. Or something like that?
    expect_true(all(abs(gsea[,'disc','avgCor',]) < .02, na.rm = TRUE))
    })

test_that('model-based singletons agree with model', {
    gsea <- gseaAfterBoot(zf, boots, sets, CoefficientHypothesis('Stim.ConditionUnstim'), control=list(n_randomize=20, var_estimate='modelbased'))@tests
    expect_equal(gsea['C','cont','var','test'], vcov(zf, 'C')['Stim.ConditionUnstim','Stim.ConditionUnstim', 15])
    expect_equal(gsea['C','disc','var','test'], vcov(zf, 'D')['Stim.ConditionUnstim', 'Stim.ConditionUnstim', 15])
})


test_that('Order is invariant', {
    setsRev <- sets[rev(seq_along(sets))]
    gseaRev <- gseaAfterBoot(zf, bootsUncor, setsRev, CoefficientHypothesis('Stim.ConditionUnstim'), control=list(n_randomize=Inf))@tests
    expect_equivalent(gsea, gseaRev[names(sets),,,])
})

gsea1 <-  gseaAfterBoot(zf, boots, sets, CoefficientHypothesis('Stim.ConditionUnstim'), control=list(n_randomize=Inf))@tests
test_that('Null and test statistic are complements in complementary modules',{
    expect_true(all.equal(gsea1['E',,c('stat', 'var', 'dof'),'null'],     gsea1['notE',,c('stat', 'var', 'dof'),'test']))
})


test_that('Variances are positive', {
     expect_true(all(gsea1[,,'var',]>0, na.rm=TRUE))
})

test_that('Avg cor is between 0 and 1 ', {
    expect_true(all(abs(gsea1[,,'avgCor',]) <= 1, na.rm = TRUE))
    })

test_that('combining coefficients works',{
              Zn <- calcZ(gseaClass, testType='normal')
              Zt <- calcZ(gseaClass, testType='t')
              ntest <- rowSums(!is.na(Zn[,,'Z']))
              ZnF <- calcZ(gseaClass, testType='normal', combined='fisher')
              ZtF <- calcZ(gseaClass, testType='t', combined='fisher')
              ZnS <- calcZ(gseaClass, testType='normal', combined='sto')
              ZtS <- calcZ(gseaClass, testType='t', combined='sto')
              ## Stouffer is just scale sum of Z's under normality
              expect_equal(ZnS[,'Z'], rowSums(Zn[,,'Z'],na.rm=TRUE)/sqrt(ntest))
              ## Fisher is just chi-square sum of log-pvalues
              expect_equal(ZnF[,'P'], pchisq(rowSums(-2*log(Zn[,,'P']),na.rm=TRUE), df=2*ntest, lower.tail=FALSE))
              expect_equal(ZtF[,'P'], pchisq(rowSums(-2*log(Zt[,,'P']),na.rm=TRUE), df=2*ntest, lower.tail=FALSE))
              ## t is strictly smaller than normal
              expect_true(all(ZnS[,'P']<ZtS[,'P']))
              expect_true(all(ZnF[,'P']<ZtF[,'P']))
              ## Stouffer T is extremely similar to normal
              expect_gt(cor(ZnS[,'Z'], ZtS[,'Z']), .99)
              ## Fisher Z scores are reasonable
              expect_gt(cor(abs(ZnF[,'Z']), -log(ZnF[,'P'])), .85)
          })

test_that('Summary method works', {
    s <- summary(gseaClass)
    })
