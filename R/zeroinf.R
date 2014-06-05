methodDict <- c('glm'='GLMlike', 'lmer'='LMERlike')

##' Convenience function for running a zero-inflated regression
##'
##' Fits a hurdle model on zero-inflated continuous data in which the zero process
##' is modeled as a logistic regression
##' and (conditional on the the response being >0), the continuous process is Gaussian, ie, a linear regression.
##' @param formula model formula
##' @param data a data.frame, list or environment in which formula is evaluated
##' @param method one of 'glm' or 'lmer'.  See SingleCellAssay:::methodDict for other possibilities.
##' @param silent if TRUE suppress common errors from fitting continuous part
##' @param ... passed to \code{fit}, and eventually to the linear model fitting function
##' @return list with "disc"rete part and "cont"inuous part 
##' @export
##' @examples
##' data<- data.frame(x=rnorm(500), z=rbinom(500, 1, .3))
##' logit.y <- with(data, x*2 + z*2); mu.y <- with(data, 10+10*x+10*z + rnorm(500))
##' y <- (runif(500)<exp(logit.y)/(1+exp(logit.y)))*1
##' y[y>0] <- mu.y[y>0]
##' data$y <- y
##' fit <- zlm(y ~ x+z, data)
##' summary(fit$disc)
##' summary(fit$cont)
##'
##' @seealso GLMlike, LMERlike
zlm <- function(formula, data, method='glm',silent=TRUE, ...){
  if(!inherits(data, 'data.frame')) stop("'data' must be data.frame, not matrix or array")
  if(!is(formula, 'formula')) stop("'formula' must be class 'formula'")

  ## lm initially just to get response vector
  ## Turn glmer grouping "|" into "+" to get correct model frame
  resp <- eval(formula[[2]], data)
  obj <- new(methodDict[method], formula=formula, design=data, response=resp)
  obj <- fit(obj)
  list(cont=obj@fitC, disc=obj@fitD)
}

summary.zlm <- function(out){
  summary(out$cont)
  summary(out$disc)
}

##' zero-inflated regression for SingleCellAssay 
##'
##' For each gene in sca, fits the hurdle model in \code{formula} (linear for et>0), logistic for et==0 vs et>0.
##' Conducts tests using hypothesis.matrix.
##'
##' When keep.zlm is FALSE, a 3D array with first dimension being the genes,
##' next dimension giving information about the test
##' (the degrees of freedom, Chisq statistic, and P value), and final dimension
##' being the value of these quantities on the
##' discrete, continuous and hurdle (combined) levels.
##'
##' When keep.zlm is TRUE, a list of length two is returned.
##' Component "tests" gives the above 3-D array.
##' Component "models" is a list giving the model fit for each gene.
##' @title zlm.SingleCellAssay
##' @param formula a formula with the measurement variable on the LHS and predictors present in cData on the RHS
##' @param sca SingleCellAssay object
##' @param method character vector, either 'glm' or 'glmer'
##' @param hypo.terms character vector giving terms to drop from model
##' @param hypo.contrasts specific contrasts to test in form expected by lht
##' @param type type of test to run, one of 'Wald' or 'LRT'
##' @param keep.zlm should the model objects be kept?
##' @param .parallel run fits using parallel processing.  must have doParallel
##' @param silent Silence common problems with fitting some genes
##' @param ... passed to glm/glmer
##' @return either an array of tests (one per primer) or a list
##' @export
##' @importFrom stringr str_split_fixed
##' @importFrom stringr fixed
##' @examples
##' \dontrun{
##' data(vbetaFA)
##' testsByGene <- zlm.SingleCellAssay(Et ~ Stim.Condition, vbetaFA, hypothesis.matrix='Stim.ConditionUnstim')
##' # genes X metric X test type
##' dimnames(testsByGene)
##'
##' modelsAndTestsByGene <- zlm.SingleCellAssay(Et ~ Stim.Condition, vbeta.sc, hypothesis.matrix='Stim.ConditionUnstim', keep.zlm=TRUE)
##' names(modelsAndTestsByGene$models)
##' summary(modelsAndTestsByGene$models[['IL13']]$disc)
##' summary(modelsAndTestsByGene$models[['IL13']]$cont)
##' }
zlm.SingleCellAssay <- function(formula, sca, method='glm', hypothesis, type='Wald', keep.zlm='false', .parallel=FALSE, silent=TRUE, ...){
    method <- match.arg(method, c('glm', 'glmer'))
    type <- match.arg(type, c('LRT', 'Wald'))
    test <- if(type=='LRT') lrTest else waldTest
    
    if(!is(sca, 'SingleCellAssay')) stop("'sca' must be (or inherit) 'SingleCellAssay'")
    if(!is(formula, 'formula')) stop("'formula' must be class 'formula'")
    fsplit <- str_split_fixed(deparse(formula), fixed('~'), 2)
    if(nchar(fsplit[1,1])>0) message("Ignoring LHS of formula (", fsplit[1,1], ') and using exprs(sca)')
    formula <- as.formula(paste0('~', fsplit[1,2]))
    obj <- new(methodDict[method], design=cData(sca), formula=formula)
    
    genes <- colnames(exprs(sca))
    ng <- length(genes)
    upperQgene <- which(rank(freq(sca))==floor(.75*ng))
    obj <- fit(obj, exprs(sca)[,upperQgene])

    testNames <- makeChiSqTable(c(0, 0), c(1, 1), '')
    coefNames <- names(coef(obj, 'C'))
    vcovNames <- colnames(vcov(obj, 'C'))
    tests <- array(0, dim=c(ng, nrow(testNames), ncol(testNames)), dimnames=list(primerid=genes, test.type=row.names(testNames), metric=colnames(testNames)))
    ## Todo: coefs, vcov, etc
    ## coef <- 
    
    for(i in seq_len(ng)){
        tt <- try({
            obj <- fit(obj, response=exprs(sca)[,i])
            test(obj, hypothesis)
        }, silent=silent)
        if(is(tt, 'try-error')) next
        tests[i,,] <- tt
    }

    structure(tests, obj=obj)
}
