
zlm <- function(formula, data, silent=TRUE, ...){
  init <- lm(formula, data, method='model.frame', ...) #run lm initially just to get pull response vector
  data$pos <-   model.response(init)>0
  cont <- try(glm(formula, data, subset=pos, family='gaussian', ...), silent=silent)
  if(inherits(cont, 'try-error')){
    warning('Some factors were not present among the positive part')
    cont <- lm(0~0)
  }                                     
  
  init[[1L]] <- (init[[1L]]>0)*1            #apparently the response goes first in the model.frame
  disc <- glm(formula, init, family='binomial')
  out <- list(cont=cont, disc=disc)
  class(out) <- 'zlm'
  out
}

summary.zlm <- function(out){
  summary(out$cont)
  summary(out$disc)
}

## Do LR test separately on continuous and discrete portions
test.zlm <- function(model, scope){
  tt <- try({
    if(missing(scope)){
      d.c <- drop1(model$cont, scope=scope, test='LRT')
} else{
 d.c <- drop1(model$cont, test='LRT')
}
  }, TRUE)
  
  if(inherits(tt, 'try-error')){
    d.c <- 0
  }

  tt <- try({
    if(missing(scope)){
      d.d <- drop1(model$disc, scope=scope, test='LRT')
} else{
 d.d <- drop1(model$disc, test='LRT')
}
  }, TRUE)
  
  if(inherits(tt, 'try-error')){
    d.d <- 0
  }

  
  
  d.n <- d.d+d.c
  
  d.n[,'Pr(>Chi)'] <- NA
  d.n
}

## terms: output from terms(formula)
## var: character of a variable that appeared in a formula (including interactions)
## mm: model matrix (output of model.matrix(formula, data))
## returns names
coefsForVar <- function(terms, mm, term){
 assignIdx <- which(labels(terms) == term) #gives us index into assign attribute
          varCoefIdx <- which(attr(mm, 'assign') == assignIdx)     #gives coefficient idx corresponding to term in formula
          coefForThisVar <- colnames(mm)[varCoefIdx]
 coefForThisVar

}

naToZero <- function(numeric){
  numeric[is.na(numeric)] <- 0
  numeric
}

##' zero-inflated regression for SingleCellAssay 
##'
##' Fits a hurdle model in \code{formula} (linear for et>0), logistic for et==0 vs et>0.
##' Conducts likelihood ratio tests for each predictor in \code{formula} that does not appear in
##' \code{scope}.
##'
##' A \code{list} of \code{data.frame}s, is returned, with one \code{data.frame} per tested predictor.
##' Rows of each \code{data.frame} are genes, the columns give the value of the LR test and its P-value, and the sum of the T-statistics for each level of the factor (when the predictor is categorical).
##' @title zlm.SingleCellAssay
##' @param formula a formula with the measurement variable on the LHS and predictors present in cData on the RHS
##' @param sca SingleCellAssay object
##' @param scope (optional) a formula giving the size of the smaller model to be fit.  If omitted, each predictor will be dropped in turn.
##' @param ... passed to lm and glm. 
##' @return a \code{list} of \code{data.frame}, one per tested predictor.  See details.
##' @export
##' @importFrom reshape rename
##' @importFrom stringr str_match
zlm.SingleCellAssay <- function(formula, sca, scope, ...){
    probeid <- getMapping(sca@mapping,"primerid")[[1]]

    m <- SingleCellAssay:::melt(sca)
    ssca <- split(m, m[,probeid], drop=TRUE)

    for( i in seq_along(ssca)){
      x <- ssca[[i]]
      this.fit <- zlm(formula, x, ...)
      z.disc <- summary(this.fit$disc)$coef[,3]
      z.cont <- try(summary(this.fit$cont)$coef[,3], silent=TRUE)

      if(inherits(z.cont, 'try-error'))
        z.cont <- rep(NA, length(z.disc))
      
      raw <- test.zlm(this.fit, scope)
      out <- raw[-1,][,c('Df', 'LRT')]
      out <- rename(out, c('LRT'='lrstat'))
      out$p.value <- pchisq(out$lrstat, df=out$Df, lower.tail=FALSE)
      coef.types <- c('disc', 'cont', 'sum') #what do we report about the regressions
      if( i == 1){                      #this will fail if the first gene was exceptional in some way, but trying to guess the dimension of the result is hard...
         mm <- model.matrix(formula, x, ...)
         message('Coefficients in model: \n', paste(colnames(mm), collapse = ' '))
        tests <- lapply(row.names(out), function(term){
         coefForThisVar <- coefsForVar(terms(this.fit$disc), mm, term)
          df.skeleton <- data.frame(c(out[1,], rep(NA, length(coef.types)*length(coefForThisVar))))
          names(df.skeleton)[-1:-length(out[1,])] <- outer(coefForThisVar, coef.types, FUN=paste, sep='.')
          cbind(primerid=names(ssca), df.skeleton[rep(1,length(ssca)),])        
      })
         names(tests) <- row.names(out)
      }
      
      for(term in row.names(out)){
        tests[[term]][i,][names(out[term,])] <- out[term,]
        coefForThisVar <- coefsForVar(terms(this.fit$disc), mm, term)
        tests[[term]][i,][paste(coefForThisVar, coef.types[1], sep='.')] <- z.disc[coefForThisVar]
        tests[[term]][i,][paste(coefForThisVar, coef.types[2], sep='.')] <- z.cont[coefForThisVar]
        tests[[term]][i,][paste(coefForThisVar, coef.types[3], sep='.')] <- naToZero(abs(z.disc[coefForThisVar])) + naToZero(abs(z.cont[coefForThisVar]))
      }
                    
    }
    
    tests
}
