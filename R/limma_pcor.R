#' Test partial correlation of each row of an object to a phenotype vector
#' 
#' Test partial correlation of each row of an object to a phenotype vector given covariates.
#' 
#' @param object A matrix-like data object containing log-ratios or 
#'   log-expression values, with rows corresponding to 
#'   genes and columns to samples.
#' @param phenotype Numeric vector of phenotypes of the samples. Should be same length as
#'   \code{ncol(object)}. If the vector is named, names should match 
#'   \code{colnames(object)}.
#' @param covariates Numeric matrix or dataframe with one or more covariates as columns passed to 
#' \code{\link[limma]{removeBatchEffect}}.
#' @param reorder.rows logical, should rows be reordered by p-value or be left in 
#'   the same order as \code{object}?
#' @param prefix character string to add to beginning of column names.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' @return Dataframe.

limma_pcor <- function(object, phenotype, covariates, reorder.rows=TRUE, prefix=NULL, adjust.method='BH'){
  stopifnot(length(phenotype)==ncol(object), names(phenotype)==colnames(object), limma::isNumeric(covariates),
            nrow(as.matrix(covariates))==ncol(object))
  
  #pheno residuals
  #data.frame(cbind()) can handle NULLs but not factors
  dat <- data.frame(cbind(phenotype, covariates))
  
  pheno.fm <- stats::lm(phenotype ~ ., data=dat)
  pheno.res <- stats::residuals(pheno.fm)
  #names get corrupted
  names(pheno.res) <- names(phenotype)
  
  #object residuals
  object.res <- limma::removeBatchEffect(x=object, covariates = covariates)
  
  #how many df to remove in limma?
  reduce.df <- ncol(as.matrix(covariates))
  
  lc <- limma_cor(object=object.res, phenotype=pheno.res, reduce.df=reduce.df, reorder.rows=reorder.rows, prefix=prefix, 
                  adjust.method=adjust.method)
  lc <- lc[,setdiff(colnames(lc), grep("\\.AveExpr$", colnames(lc), value=TRUE))]
  return(lc)
}
