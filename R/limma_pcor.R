#' Test partial correlation of each row of an object to a phenotype vector
#' 
#' Test partial correlation of each row of an object to a phenotype vector given covariates.
#' 
#' @inheritParams limma_cor
#' @inheritParams hitman
#' @inheritParams limma_contrasts
#' @return Data frame.

limma_pcor <- function(object, phenotype, covariates, reorder.rows=TRUE, prefix=NULL, adjust.method="BH", 
                       check.names=TRUE, cols=c("t", "P.Value")){
  stopifnot(length(phenotype)==ncol(object), limma::isNumeric(covariates), nrow(as.matrix(covariates))==ncol(object))
  if (check.names){ stopifnot(names(phenotype)==colnames(object)) }
  
  # pheno residuals
  # data.frame(cbind()) can handle NULLs but not factors
  dat <- data.frame(cbind(phenotype, covariates))
  pheno.fm <- stats::lm(phenotype ~ ., data=dat)
  pheno.res <- stats::residuals(pheno.fm)
  # names get corrupted
  names(pheno.res) <- names(phenotype)

  # object residuals
  object.res <- limma::removeBatchEffect(x=object, covariates = covariates)
  
  # how many df to remove in limma?
  reduce.df <- ncol(as.matrix(covariates))
  lc <- limma_cor(object=object.res, phenotype=pheno.res, reduce.df=reduce.df, reorder.rows=reorder.rows, prefix=prefix, 
                  adjust.method=adjust.method, cols=cols)
  return(lc)
}
