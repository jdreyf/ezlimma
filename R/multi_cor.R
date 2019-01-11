#' Test correlation of each row of an object to each column of pheno.mat
#' 
#' Test correlation of each row of an object to each column of \code{pheno.mat} using 
#' one of Pearson's, Kendall's, or Spearman's correlation methods, or limma 
#' regression in \code{\link{limma_cor}}.
#' 
#' @param object A matrix-like data object containing log-ratios or 
#'   log-expression values, with rows corresponding to features (e.g. genes) and 
#'   columns to samples.
#' @param pheno.mat matrix-like data object of phenotypes of the samples, with each column one 
#'   phenotype vector. Length and names of rows of \code{pheno.mat} should 
#'   correspond to columns of \code{object}.
#' @param method a character string indicating which association is to be used 
#'   for the test. One of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, 
#'   from \code{\link[stats]{cor.test}} or \code{"limma"} for \code{\link{limma_cor}}.
#' @param reorder.rows logical, should rows be reordered by F-statistic from 
#'   \code{\link[limma]{toptable}} or be left in the same order as 
#'   \code{object}?
#' @param prefix character string to add to beginning of column names.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' @param limma.cols if \code{method="limma"}, this specifies \code{cols} from 
#'   \code{\link{limma_cor}}. Ignored without a warning if \code{method} not 
#'   \code{"limma"}.
#' @return Dataframe with several statistical columns corresponding to each
#'   phenotype and one row per feature.
#' @details  Each column of \code{pheno.mat} is tested independently.
#' @export

multi_cor <- function(object, pheno.mat, method=c('pearson', 'spearman', 'kendall', 'limma'),
                     reorder.rows=TRUE, prefix=NULL, adjust.method='BH', 
                     limma.cols=c('AveExpr', 'P.Value', 'adj.P.Val', 'logFC')){
  method <- match.arg(method)
  if (is.null(dim(pheno.mat))) stop("pheno.mat needs to have rows and columns.")
  stopifnot(ncol(object)==nrow(pheno.mat), rownames(pheno.mat)==colnames(object))
  cor.mat <- NULL
  for (i in 1:ncol(pheno.mat)){
    prefix.tmp <- ifelse(!is.null(prefix), paste(prefix, colnames(pheno.mat)[i], sep='.'), colnames(pheno.mat)[i])
    if (method=='limma'){
      cor.tmp <- data.matrix(limma_cor(object, pheno.mat[,i], reorder.rows=FALSE, prefix=prefix.tmp, cols=limma.cols))
    } else {
      cor.tmp <- ezcor(object, pheno.mat[,i], method=method, reorder.rows=FALSE, prefix=prefix.tmp, adjust.method=adjust.method)
    }
    if (!is.null(cor.mat)){ stopifnot(rownames(cor.mat)==rownames(cor.tmp)) }
    cor.mat <- cbind(cor.mat, cor.tmp)
  }
  if (reorder.rows){ cor.mat <- cor.mat[order(combine_pvalues(cor.mat)),] }
  return(cor.mat)
}