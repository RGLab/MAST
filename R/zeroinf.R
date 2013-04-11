
zlm <- function(formula, data, ...){
  init <- lm(formula, data, method='model.frame', ...) #run lm initially just to get pull response vector
  data$pos <-   model.response(init)>0
  cont <- try(glm(formula, data, subset=pos, family='gaussian', ...))
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

zlm.SingleCellAssay <- function(formula, sca, scope, ...){
    probeid <- getMapping(sca@mapping,"primerid")[[1]]
    measure <- getMapping(sca@mapping,"measurement")[[1]]

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
      if( i == 1){                      #this will fail if the first gene was exceptional in some way, but trying to guess the dimension of the result is hard...
        tests <- lapply(1:nrow(out), function(vIdx){
          var <- labels(terms(this.fit$disc))[vIdx]
          coefNames <- names(coef(this.fit$disc))
          coefMatch <- str_match(pattern=var, string=coefNames) # covar.names x coefficient lookup table
          coefForThisVar <- coefNames[!is.na(coefMatch)]
          df.skeleton <- data.frame(c(out[1,], rep(NA, length(coefForThisVar))))
          names(df.skeleton)[-1:-length(out[1,])] <- coefForThisVar
          cbind(gene=names(ssca), df.skeleton[rep(1,length(ssca)),])        
      })
        names(tests) <- labels(terms(this.fit$disc))
      }
      
      for(vIdx in seq_len(nrow(out))){
        tests[[vIdx]][i,][names(out[vIdx,])] <- out[vIdx,]
        var <- labels(terms(this.fit$disc))[vIdx]
         coefNames <- names(coef(this.fit$disc))
          coefMatch <- str_match(pattern=var, string=coefNames) # covar.names x coefficient lookup table
          coefForThisVar <- coefNames[!is.na(coefMatch)]
        tests[[vIdx]][i,][coefForThisVar] <- z.disc[coefForThisVar]+z.cont[coefForThisVar]
      }
                    
    }
    
    tests
}
