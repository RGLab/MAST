#' @name CovFromBoots
#' @title Extract the inter-gene covariance matrices for continuous and discrete components of a MAST model for a given coefficient from bootstrap replicates
#' @description Computes the genewise covariance for a model coefficient from bootstrap replicates
#' from `MAST::bootVcov1()`. If coefficients are unestimable (i.e. NA) for a gene, that row/column in the
#' covariance matrix will be NA. Returns a list with components "C" and "D" containing the covariance
#' matrices for the "C"ontinuous and "D"iscrete components of the MAST model.
#' @return list with components "C" and "D" containing covariance matrices for the continuous and discrete components of the model.
#' @param boots a multidimensional array returned by `bootVcov1` or `pbootVcov1`.
#' @param coefficient `character` the name of the model coefficient for which to return the inter-gene covariance matrices.
CovFromBoots <- function(boots = NULL,
								 coefficient = NULL) {
	if (is.null(coefficient)) {
		stop("You must specify a model coefficient.")
	}
	if (is.null(boots)) {
		stop("You must provide bootstrap results.")
	}
	pboots <- aperm(boots, perm = c(2, 3, 4, 1))
	bmean <- rowMeans(pboots, dims = 3)
	l <- apply(pboots, 4, function(x) {
		xc <- x - bmean
		dim(xc) <- dim(x)
		list(xc)
	})
	l <- lapply(l, function(x)
		x[[1]])
	cboots <- abind::abind(l, along = 4)
	
	dimnames(cboots)[1:3] <- dimnames(bmean)
	dimnames(cboots) <-
		setNames(dimnames(cboots), c('genes', 'coef', 'comp', 'rep'))
	
	
	bootstat <- cboots[, coefficient, , ]
	
	naboots <- rowSums(is.na(bootstat), dims = 2)
	naAny <- naboots
	
	
	bootstat[is.na(bootstat)] <- 0
	
	covariances <- list(C = NULL, D = NULL)
	for (component in c("C", "D")) {
		i <- bootstat[, , , drop = FALSE]
		j <- bootstat[, , , drop = FALSE]
		
		tci <-
			(function(component = component,
						 i = i) {
				sub <-
					i[, component, , drop = FALSE]
				dim(sub) <- dim(sub)[-2]
				sub
			})(component, i)
		tcj <-
			(function(component = component,
						 j = j) {
				sub <-
					j[, component, , drop = FALSE]
				dim(sub) <- dim(sub)[-2]
				sub
			})(component, j)
		tcp <- tcrossprod(tci, tcj)
		tcp[(naAny > 0)[, component], ] <- NA
		tcp[, (naAny > 0)[, component]] <- NA
		covariances[[component]] <- tcp / (dim(tci)[length(dim(tci))] - 1)
	}
	covariances
}
