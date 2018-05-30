#' Test partial correlation of each row of an object to a phenotype vector
#' 
#' Test partial correlation of each row of an object to a phenotype vector given covariates.
#' 
#' @param object A matrix-like data object containing log-ratios or 
#'   log-expression values, with rows corresponding to 
#'   genes and columns to samples.
#' @param phenotype Vector of phenotypes of the samples. Should be same length as
#'   \code{ncol(object)}. If the vector is named, names should match 
#'   \code{colnames(object)}.
#' @param batch Batch vector passed to \code{\link[limma]{removeBatchEffect}}, which coerces it to factor.
#' @param covariates Numeric matrix with one or more covariates as columns passed to 
#' \code{\link[limma]{removeBatchEffect}}.
#' @param reorder.rows logical, should rows be reordered by p-value or be left in 
#'   the same order as \code{object}?
#' @param prefix character string to add to beginning of column names.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' @return Dataframe.
#' @details \code{\link[limma]{removeBatchEffect}} treats \code{batch} differently form \code{covariates}, even if 
#' \code{covariates=as.numeric(as.factor(batch))}.
#' @import stats
#' @importFrom limma removeBatchEffect

limma_pcor <- function(object, phenotype, batch=NULL, covariates=NULL, reorder.rows=TRUE, prefix=NULL, adjust.method='BH'){
  stopifnot(length(phenotype)==ncol(object), names(phenotype)==colnames(object), !is.null(batch)|!is.null(covariates),
            is.null(batch)|length(batch)==ncol(object),
            is.null(covariates)|is.numeric(covariates))
  
  #pheno residuals
  if (is.null(covariates)){
    pheno.fm <- stats::lm(phenotype ~ ., data=data.frame(phenotype, batch))
  } 
  if (is.null(batch)){
    pheno.fm <- stats::lm(phenotype ~ ., data=data.frame(phenotype, covariates))
  }
  if (!is.null(batch) & !is.null(covariates)){
    pheno.fm <- stats::lm(phenotype ~ ., data=data.frame(phenotype, batch, covariates))
  }
  pheno.res <- stats::residuals(pheno.fm)
  #names get corrupted
  names(pheno.res) <- names(phenotype)
  
  #object residuals
  object.res <- limma::removeBatchEffect(x=object, batch=batch, covariates = covariates)
  
  #how many df to remove in limma?
  reduce.df <- ifelse(test=is.null(covariates), yes=1, no=1+ncol(as.matrix(covariates)))
  
  lc <- limma_cor(object=object.res, phenotype=pheno.res, reduce.df=reduce.df, reorder.rows=reorder.rows, prefix=prefix, 
                  adjust.method=adjust.method)
  return(lc)
}
