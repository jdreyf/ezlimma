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
#' @param covar A matrix or dataframe with one or more covariates as columns.
#' @param reorder.rows logical, should rows be reordered by p-value or be left in 
#'   the same order as \code{object}?
#' @param prefix character string to add to beginning of column names.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' @return Dataframe.
#' @export
#' @import stats
#' @importFrom limma removeBatchEffect

ezpcor <- function(object, phenotype, covar, reorder.rows=TRUE, prefix=NULL, adjust.method='BH'){
  stopifnot(length(phenotype)==ncol(object), names(phenotype)==colnames(object), nrow(covar)==ncol(object))
  
  #pheno residuals
  pheno.fm <- stats::lm(phenotype ~ ., data=data.frame(phenotype, covar))
  pheno.res <- stats::residuals(pheno.fm)
  
  #object residuals
  object.res <- limma::removeBatchEffect(x=object, covariates = covar)
  
  #how many df to remove in limma?
  reduce.df <- ncol(as.matrix(covar))
  
  lc <- limma_cor(object=object.res, phenotype=pheno.res, reduce.df=reduce.df, reorder.rows=reorder.rows, prefix=prefix, 
                  adjust.method=adjust.method)
  return(lc)
}
